function T = hosvd(X,tol,varargin)
%HOSVD Compute sequentially-truncated higher-order SVD (Tucker).
%
%   T = HOSVD(X,TOL) computes a Tucker decomposition with relative error
%   specified by TOL, i.e., it computes a ttensor T such that 
%
%              ||X-T||/||X|| <= TOL. 
% 
%   The method automatically determines the appropriate ranks of the
%   Tucker decomposition. By default, the method computes the
%   sequentially-truncated HOSVD. 
%
%   T = HOSVD(X,TOL,'param',value,...) specifies optional parameters and
%   values. Valid parameters and their default values are:
%      'verbosity' - How much to print between 0 and 10. Default: 1.
%      'dimorder' - Order to loop through dimensions Default: 1:ndims(X).
%      'sequential' - Use sequentially-truncated version: Default: true.
%      'ranks' - Specify ranks (rather than computing). Default: [].
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html',hosvd_doc.html')))">Documentation page for HOSVD</a>
%
%   See also TUCKER_ALS, TTENSOR
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

%% Read paramters
d = ndims(X);

params = inputParser;
params.addParameter('verbosity',1);
params.addParameter('sequential',true);
params.addParameter('dimorder',1:d,@(x) isequal(sort(x),1:d));
params.addParameter('ranks',[]);
params.parse(varargin{:});

verbosity = params.Results.verbosity;
sequential = params.Results.sequential;
dimorder = params.Results.dimorder;
ranks = params.Results.ranks;

%% Setup
if verbosity > 0
    fprintf('Computing HOSVD...\n');
end
normxsqr = collapse(X.^2);
eigsumthresh = tol.^2 * normxsqr / d;

if verbosity > 2
    fprintf('||X||^2 = %g\n', normxsqr);
    fprintf('tol = %g\n', tol);
    fprintf('eigenvalue sum threshold = tol^2 ||X||^2 / d = %g\n', eigsumthresh);
end

if ~isempty(ranks)
    if ~isvector(ranks) || length(ranks) ~= d
        error('Specified ranks must be a vector of length ndims(X)');
    end
    r = ranks;
else
    r = zeros(d,1);
end

%% Main loop

U = cell(d,1); % Allocate space for factor matrices
Y = X; % Copy input tensor, shrinks at each step for sequentially-truncated

for k = dimorder
    
    % Compute Gram matrix
    Yk = double(tenmat(Y,k));
    Z = Yk*Yk';
    
    % Compute eigenvalue decomposition
    [V,D] = eig(Z);
    [eigvec,pi] = sort(diag(D),'descend');
    
    % If rank is not prespecified, compute it.
    if r(k) == 0
        
        eigsum = cumsum(eigvec,'reverse');
        r(k) = find(eigsum > eigsumthresh, 1, 'last');
        
        if verbosity > 5
            fprintf('Reverse cummulative sum of evals of Gram matrix:\n');
            for i = 1:length(eigsum)
                fprintf('%d: %6.4f ',i,eigsum(i));
                if i == r(k)
                    fprintf('<-- Cutoff');
                end
                fprintf('\n');
            end
        end
        
    end
    
    % Extract factor matrix by picking out leading eigenvectors of V
    U{k} = V(:,pi(1:r(k)));
    
    % Shrink!
    if sequential
        Y = ttm(Y,U{k}',k);
    end
end

% Extract final core
if sequential
    G = Y;
else
    G = ttm(Y,U,'t');
end

%% Final result
T = ttensor(G,U);

if verbosity > 0
    diffnormsqr = collapse((X-full(T)).^2);
    relnorm = sqrt(diffnormsqr/normxsqr);
    fprintf('Size of core: %s\n', tt_size2str(size(G)));
    if relnorm <= tol
        fprintf('||X-T||/||X|| = %g <=', relnorm);
        fprintf('%f (tol)\n',tol);
    else
        fprintf('Tolerance not satisfied!! ');
        fprintf('||X-T||/||X|| = %g >=', relnorm);
        fprintf('%f (tol)\n',tol);
        warning('Specified tolerance was not achieved');
    end
    fprintf('\n');
end
