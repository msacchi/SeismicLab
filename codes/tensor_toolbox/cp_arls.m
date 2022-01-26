function [P,Uinit,output] = cp_arls(X,R,varargin)
%CP_ARLS CP decomposition of dense tensor via randomized least squares.
%
%   M = CP_ARLS(X,R) computes an estimate of the best rank-R
%   CP model of a dense tensor X using a randomized alternating
%   least-squares algorithm. The input X must be a (dense) tensor. The
%   result P is a ktensor. 
%
%   *Important Note:* The fit computed by CP_ARLS is an approximate
%   fit, so the stopping conditions are necessarily more conservative. The
%   approximation is based on sampling entries from the full tensor and
%   estimating the overall fit based on their individual errors.
%
%   M = CP_ARLS(X,R,'mix',0) skips the 'mixing' which is an expensive
%   preprocessing step.  In many cases, this step is not necessary and
%   requires less initialization time and space.  It is suggested to try
%   this out.
%
%   M = CP_ARLS(X,R,'param',value,...) specifies optional parameters and
%   values. Valid parameters and their default values are:
%   o 'mix'       - Include FJLT transformations {true}
%   o 'epoch'     - Number of iterations between convergence checks {50}
%   o 'maxepochs' - Maximum number of epochs {1000}
%   o 'newitol'   - Quit after this many epochs with no improvement {5}
%   o 'tol'       - Tolerance for improvement, i.e., fit - maxfit > tol {0}
%   o 'fitthresh' - Terminate when fit > fitthresh {1.000}
%   o 'printitn'  - Print fit every n epochs; 0 for no printing {10}
%   o 'init'      - Initial guess ['random'|'nvecs'|cell array] {random}
%   o 'nsamplsq'  - Number of least-squares row samples {10Rlog2(R)}
%   o 'nsampfit'  - Number of entry samples for approximate fit {2^14}
%   o 'dimorder'  - Order to loop through dimensions {1:ndims(A)}
%
%   [M,U0] = CP_ARLS(...) also returns the initial guess.
%
%   [M,U0,out] = CP_ARLS(...) also returns additional output that
%   contains the input parameters and other information.
%
%   Examples:
%   info = create_problem('Size',[100 100 100],'Num_Factors',2);
%   M = cp_arls(info.Data,2);
%
%   REFERENCE: C. Battaglino, G. Ballard, T. G. Kolda, A Practical
%   Randomized CP Tensor Decomposition, SIAM J. Matrix Analysis and
%   Applications, 39(2):876-901, 2018, https://doi.org/10.1137/17M1112303.
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','cp_arls_doc.html')))">Documentation page for CP-ARLS</a>
%
%   See also CP_ALS, KTENSOR, TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%% Extract some sizes, etc.
N = ndims(X);
sz = size(X);
num_elements = prod(sz);

%% Parse parameters
params = inputParser;
params.addParameter('init', 'random', @(x) (iscell(x) || ismember(x,{'random'})));
params.addParameter('dimorder', 1:N, @(x) isequal(sort(x),1:N));
params.addParameter('printitn', 10, @isscalar);
params.addParameter('mix', true, @islogical);
params.addParameter('nsamplsq', ceil(10*R*log2(R)));
params.addParameter('maxepochs', 1000);
params.addParameter('nsampfit', 2^14);
params.addParameter('tol', 0, @isscalar);
params.addParameter('fitthresh', 1, @(x) isscalar(x) & x > 0 & x <= 1);
params.addParameter('epoch', 50)
params.addParameter('newitol', 5);
params.addParameter('truefit', false, @islogical);
params.parse(varargin{:});

% Copy from params object
fitchangetol = params.Results.tol;
maxepochs = params.Results.maxepochs;
dimorder = params.Results.dimorder;
init = params.Results.init;
printitn = params.Results.printitn;
fitthresh = params.Results.fitthresh;   % cprand will terminate if this fit is reached (default 1)
do_fft = params.Results.mix;
newitol = params.Results.newitol;
nsamplsq = params.Results.nsamplsq;
nsampfit = params.Results.nsampfit;
epochsize = params.Results.epoch;
truefit = params.Results.truefit;
    
%% Set up initial guess for U (factor matrices)
if iscell(init)
    Uinit = init;
    if numel(Uinit) ~= N
        error('OPTS.init does not have %d cells',N);
    end
    for n = dimorder(2:end)
        if ~isequal(size(Uinit{n}),[size(X,n) R])
            error('OPTS.init{%d} is the wrong size',n);
        end
    end
else
    % Observe that we don't need to calculate an initial guess for the
    % first index in dimorder because that will be solved for in the first
    % inner iteration.
    if strcmp(init,'random')
        Uinit = cell(N,1);
        for n = dimorder(2:end)
            Uinit{n} = rand(sz(n),R);
        end
    else
        error('The selected initialization method is not supported');
    end
end

%% Set up for iterations - initializing U and the fit.
U = Uinit;
U_mixed = Uinit;
diag_flips = [];

if printitn>0    
    if(do_fft)
        fprintf('\nCP-ARLS (with mixing): \n');
    else
        fprintf('\nCP_ARLS (without mixing): \n');
    end
end

%% Sample input tensor for stopping criterion
nsampfitles = min(num_elements,nsampfit);
Xfit_subs = sample_all_modes(nsampfitles, sz);
Xfit_vals = X(Xfit_subs);
Xvalsqr_mean = mean((Xfit_vals).^2);
normX = sqrt(Xvalsqr_mean * num_elements); % Approximate!

%% Mixing tensor (if needed)
if (do_fft) % --- with mixing ---
    % Compute random diagonal D_n for each factor
    diag_flips = cell(N,1);
    for n = 1:N
        diag_flips{n} = (rand(sz(n),1)<0.5)*2-1;
    end    
    
    % Extract dense data array for use with fft and bsx commands
    X_mixed = X.data;
    % Mixing is equivalent to a series of TTMs with D_n, F_n
    % However, we can use bsxfun and fft to avoid matricizing.
    for n = N:-1:1
        % Reshape the diagonal flips into a 1*...*sz(n)*...*1 tensor
        % This lets us use bsxfun along the nth dimension.
        bsxdims = ones(1,N);
        bsxdims(n) = sz(n);
        % Note that the next line arranges the flips array to work with the
        % tensor in its default memory layout 
        flips = reshape(diag_flips{n},bsxdims);
        % fft(...,[],n) operates fiber-wise on dimension n
        X_mixed = fft(bsxfun(@times, X_mixed, flips),[],n);
    end
    X_mixed = tensor(X_mixed);
else % --- without mixing ---
    X_mixed = X; % no mixing
end

%% Mixing factors (if needed)
if (do_fft)
    % Mix factor matrices: U{i} = F{i}*D{i}*U{i}
    for i = 2:N
        U_mixed{i} = fft(bsxfun(@times,U{i},diag_flips{i}));
    end
end

%% Main Loop: Iterate until convergence
maxfit = 0;
newi = 0; % number of epochs without improvement

% ALS Loop
for epoch = 1:maxepochs
  
    % Do a bunch of iterations within each epoch
    for eiters = 1:epochsize
        
        % Iterate over all N modes of the tensor
        for n = dimorder(1:end)
            
            mixinfo.dofft = do_fft;
            mixinfo.signflips = diag_flips;
            Unew = dense_sample_mttkrp(X_mixed,U_mixed,n,nsamplsq,mixinfo);
            
            if issparse(Unew)
                Unew = full(Unew);   % for the case R=1
            end
            
            % Normalize each vector to prevent singularities in coefmatrix
            if epoch == 1
                lambda = sqrt(sum(abs(Unew).^2,1))'; %2-norm
            else
                lambda = max( max(abs(Unew),[],1), 1 )'; %max-norm
            end
            
            Unew = bsxfun(@rdivide, Unew, lambda');
            U_mixed{n} = Unew;
            if (do_fft)
                U{n} = real(bsxfun(@times, ifft(Unew), diag_flips{n}));
            else
                U{n} = Unew;
            end
        end
    end
    
    % After each epoch, check convergence conditions
    P = ktensor(lambda, U);
    Pfit_vals = sample_ktensor(P, Xfit_subs);
    elem_mean = mean((Xfit_vals - Pfit_vals).^2);
    normDiff = sqrt(elem_mean * num_elements); % Approximate!
    fit = 1 -  normDiff / normX;
    
    if fit > maxfit + fitchangetol
        newi = 0;
        maxfit = fit;
        Psave = P; % Keep the best one seen so far!
    else
        newi = newi + 1;
    end
    
    if (fit > fitthresh) || (newi >= newitol)
        flag = 0;
    else
        flag = 1;
    end
    
    if (mod(epoch,printitn)==0) || ((printitn>0) && (flag==0))
        fprintf(' Iter %2dx%d: f~ = %e newi = %d\n', epoch, epochsize, fit, newi);
    end
    
    % Check for convergence
    if (flag == 0)
        break;
    end
end
%% Clean up final result
% Arrange the final tensor so that the columns are normalized.
P = Psave;
P = arrange(P);
P = fixsigns(P); % Fix the signs

if truefit
    normresidual = sqrt( norm(X)^2 + norm(P)^2 - 2 * innerprod(X,P) );
    fit = 1 - (normresidual / normX);%fraction explained by model
    Pfit_vals = sample_ktensor(P, Xfit_subs);
    Xfit_mean = mean((Xfit_vals - Pfit_vals).^2);
    testfit = 1 - sqrt(Xfit_mean*num_elements)/normX;
    if printitn > 0
        fprintf(' Final fit = %e Final estimated fit = %e \n', fit, testfit);
    end
else
    fit = NaN;
end

output = struct;
output.params = params.Results;
output.iters = epoch;
output.fit = fit;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Sub-functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% MTTKRP That performs sampling after transforming X and the KR-product with an FFT
% And then solves using normal equations
function [V, Xsamp, Zsamp] = dense_sample_mttkrp(X,U,n,nsamplsq,mixinfo)

N = ndims(X);

if (N < 2)
    error('MTTKRP is invalid for tensors with fewer than 2 dimensions');
end

if (length(U) ~= N)
    error('Cell array is the wrong length');
end

if n == 1
    R = size(U{2},2);
else
    R = size(U{1},2);
end

for i = 1:N
    if i == n, continue; end
    if (size(U{i},1) ~= size(X,i)) || (size(U{i},2) ~= R)
        error('Entry %d of cell array is wrong size', i);
    end
end

dims = size(X);

% Compute uniform samples for tensor and factor matrices
[tensor_idx, factor_idx] = sample_mode_n(nsamplsq, dims, n);

% Reshape the sampled tensor
Xsamp = reshape(X(tensor_idx), dims(n), []);

% Perform a sampled KRP
Zsamp = skr(U{[1:n-1,n+1:N]}, factor_idx);

% alternative
V = Xsamp / Zsamp.';
if (mixinfo.dofft)
    V = real(bsxfun(@times,ifft(V),mixinfo.signflips{n}));
    V = fft(bsxfun(@times,V,mixinfo.signflips{n}));
end

return;
end


% Sample Khatri-Rao Product of a cell array of factors
% Without forming the full KR Product
function P = skr(varargin)
if iscell(varargin{1}) % Input is a single cell array
    A = varargin{1};
else % Input is a sequence of matrices
    A = varargin(1:end-1);
end

numfactors = length(A);
matorder = numfactors:-1:1;
idxs = varargin{end};

%% Error check on matrices and compute number of rows in result
ndimsA = cellfun(@ndims, A);
if(~all(ndimsA == 2))
    error('Each argument must be a matrix');
end

ncols = cellfun(@(x) size(x, 2), A);
if(~all(ncols == ncols(1)))
    error('All matrices must have the same number of columns.');
end

P = A{matorder(1)}(idxs(:,matorder(1)),:);
for i = matorder(2:end)
    %P = P .*A{i}(idxs(:,i),:);
    P = bsxfun(@times, P, A{i}(idxs(:,i),:));
end
end


% Random sample fibers in mode n from tensor X
% Generate the corresponding indices for the factor matrices as a tuple
function [tensor_idx, factor_idx] = sample_mode_n(nsamplsq, dims, n)
D = length(dims);
tensor_idx = zeros(nsamplsq, D);     % Tuples that index fibers in original tensor

tensor_idx(:,n) = ones(nsamplsq, 1);
for i = [1:n-1,n+1:D]
    % Uniformly sample w.r. in each dimension besides n
    tensor_idx(:,i) = randi(dims(i), nsamplsq, 1);
end

% Save indices to sample from factor matrices
factor_idx = tensor_idx(:,[1:n-1,n+1:D]);

% Expand tensor_idx so that every fiber element is included
%tensor_idx = repelem(tensor_idx,dims(n),1); % not portable
tensor_idx = kron(tensor_idx,ones(dims(n),1)); % portable
tensor_idx(:,n) = repmat((1:dims(n))',nsamplsq,1);
tensor_idx = tt_sub2ind(dims, tensor_idx);
end


% Random sample fibers in mode n from tensor X
% Generate the corresponding indices for the factor matrices as a tuple
function [subs, idxs] = sample_all_modes(nsamplsq, dims)
D = length(dims);
subs = zeros(nsamplsq, D);     % Tuples that index fibers in original tensor

for i = 1:D
    % Uniformly sample w.r. in each dimension
    subs(:,i) = randi(dims(i), nsamplsq, 1);
end

% subs can be used to sample from factor matrices as well as the tensor
subs = unique(subs,'rows'); %todo: do this more efficiently
idxs = tt_sub2ind(dims, subs);
end


% Random sample fibers in mode n from tensor X
% Generate the corresponding indices for the factor matrices as a tuple
function [data] = sample_ktensor(P, subs)
data = skr(P.u, subs) * P.lambda;
end
