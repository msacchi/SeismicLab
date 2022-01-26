function [P,Uinit,output] = cp_als(X,R,varargin)
%CP_ALS Compute a CP decomposition of any type of tensor.
%
%   M = CP_ALS(X,R) computes an estimate of the best rank-R
%   CP model of a tensor X using an alternating least-squares
%   algorithm.  The input X can be a tensor, sptensor, ktensor, or
%   ttensor. The result M is a ktensor.
%
%   M = CP_ALS(X,R,'param',value,...) specifies optional parameters and
%   values. Valid parameters and their default values are:
%      'tol' - Tolerance on difference in fit {1.0e-4}
%      'maxiters' - Maximum number of iterations {50}
%      'dimorder' - Order to loop through dimensions {1:ndims(A)}
%      'init' - Initial guess [{'random'}|'nvecs'|cell array]
%      'printitn' - Print fit every n iterations; 0 for no printing {1}
%      'fixsigns' - Call fixsigns at end of iterations {true}
%
%   [M,U0] = CP_ALS(...) also returns the initial guess.
%
%   [M,U0,out] = CP_ALS(...) also returns additional output that contains
%   the input parameters.
%
%   Note: The "fit" is defined as 1 - norm(X-full(M))/norm(X) and is
%   loosely the proportion of the data described by the CP model, i.e., a
%   fit of 1 is perfect.
%
%   NOTE: Updated in various minor ways per work of Phan Anh Huy. See Anh
%   Huy Phan, Petr Tichavský, Andrzej Cichocki, On Fast Computation of
%   Gradients for CANDECOMP/PARAFAC Algorithms, arXiv:1204.1586, 2012.
%
%   Examples:
%   X = sptenrand([5 4 3], 10);
%   M = cp_als(X,2);
%   M = cp_als(X,2,'dimorder',[3 2 1]);
%   M = cp_als(X,2,'dimorder',[3 2 1],'init','nvecs');
%   U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of M
%   [M,U0,out] = cp_als(X,2,'dimorder',[3 2 1],'init',U0);
%   M = cp_als(X,2,out.params); %<-- Same params as previous run
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','cp_als_doc.html')))">Documentation page for CP-ALS</a>
%
%   See also KTENSOR, TENSOR, SPTENSOR, TTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%% Extract number of dimensions and norm of X.
N = ndims(X);
normX = norm(X);

%% Set algorithm parameters from input or by using defaults
params = inputParser;
params.addParameter('tol',1e-4,@isscalar);
params.addParameter('maxiters',50,@(x) isscalar(x) & x > 0);
params.addParameter('dimorder',1:N,@(x) isequal(sort(x),1:N));
params.addParameter('init', 'random', @(x) (iscell(x) || ismember(x,{'random','nvecs'})));
params.addParameter('printitn',1,@isscalar);
params.addParameter('fixsigns',true,@islogical);
params.parse(varargin{:});

%% Copy from params object
fitchangetol = params.Results.tol;
maxiters = params.Results.maxiters;
dimorder = params.Results.dimorder;
init = params.Results.init;
printitn = params.Results.printitn;

%% Error checking 

%% Set up and error checking on initial guess for U.
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
            Uinit{n} = rand(size(X,n),R);
        end
    elseif strcmp(init,'nvecs') || strcmp(init,'eigs') 
        Uinit = cell(N,1);
        for n = dimorder(2:end)
            Uinit{n} = nvecs(X,n,R);
        end
    else
        error('The selected initialization method is not supported');
    end
end

%% Set up for iterations - initializing U and the fit.
U = Uinit;
fit = 0;

% Store the last MTTKRP result to accelerate fitness computation.
U_mttkrp = zeros(size(X, dimorder(end)), R);

if printitn>0
  fprintf('\nCP_ALS:\n');
end

%% Main Loop: Iterate until convergence

if (isa(X,'sptensor') || isa(X,'tensor')) && (exist('cpals_core','file') == 3)
 
    %fprintf('Using C++ code\n');
    [lambda,U] = cpals_core(X, Uinit, fitchangetol, maxiters, dimorder);
    P = ktensor(lambda,U);
    
else
    
    UtU = zeros(R,R,N);
    for n = 1:N
        if ~isempty(U{n})
            UtU(:,:,n) = U{n}'*U{n};
        end
    end
    
    for iter = 1:maxiters
        
        fitold = fit;
        
        % Iterate over all N modes of the tensor
        for n = dimorder(1:end)
            
            % Calculate Unew = X_(n) * khatrirao(all U except n, 'r').
            Unew = mttkrp(X,U,n);
            % Save the last MTTKRP result for fitness check.
            if n == dimorder(end)
              U_mttkrp = Unew;
            end
            
            % Compute the matrix of coefficients for linear system
            Y = prod(UtU(:,:,[1:n-1 n+1:N]),3);
            Unew = Unew / Y;
            if issparse(Unew)
                Unew = full(Unew);   % for the case R=1
            end
                        
            % Normalize each vector to prevent singularities in coefmatrix
            if iter == 1
                lambda = sqrt(sum(Unew.^2,1))'; %2-norm
            else
                lambda = max( max(abs(Unew),[],1), 1 )'; %max-norm
            end            
            
            Unew = bsxfun(@rdivide, Unew, lambda');

            U{n} = Unew;
            UtU(:,:,n) = U{n}'*U{n};
        end
        
        P = ktensor(lambda,U);

        % This is equivalent to innerprod(X,P).
        iprod = sum(sum(P.U{dimorder(end)} .* U_mttkrp) .* lambda');
        if normX == 0
            fit = norm(P)^2 - 2 * iprod;
        else
            normresidual = sqrt( normX^2 + norm(P)^2 - 2 * iprod );
            fit = 1 - (normresidual / normX); %fraction explained by model
        end
        fitchange = abs(fitold - fit);
        
        % Check for convergence
        if (iter > 1) && (fitchange < fitchangetol)
            flag = 0;
        else
            flag = 1;
        end
        
        if (mod(iter,printitn)==0) || ((printitn>0) && (flag==0))
            fprintf(' Iter %2d: f = %e f-delta = %7.1e\n', iter, fit, fitchange);
        end
        
        % Check for convergence
        if (flag == 0)
            break;
        end        
    end   
end


%% Clean up final result
% Arrange the final tensor so that the columns are normalized.
P = arrange(P);
% Fix the signs
if params.Results.fixsigns
    P = fixsigns(P);
end

if printitn>0
    if normX == 0
        fit = norm(P)^2 - 2 * innerprod(X,P);
    else
        normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
        fit = 1 - (normresidual / normX); %fraction explained by model
    end
  fprintf(' Final f = %e \n', fit);
end

output = struct;
output.params = params.Results;
output.iters = iter;
