function data = fg_setup(Model, A, varargin)
%FG_SETUP Setup for optimization of symmetric Kruskal model.
%
%   DATA = FG_SETUP(MODEL,A,'param',value,...) is a setup routine
%   for computing an objective function that compares a symmetric
%   Kruskal tensor, MODEL, with a symmetric tensor, A. Parameter-value
%   pairs control the exact definition of the objective function and
%   constraints. For an M-way, N-dimension tensor A, a symmetric Kruskal
%   tensor model is defined by a P-vector LAMBDA and an NxP matrix X, where
%   P is the number of components. The number of optimization variables is
%   Q = P*(N+1). The DATA produces by this setup function can be used
%   repeatedly for the same tensor, A, and *any* model that is the same
%   size as MODEL,  i.e., the same number of components, same number of
%   modes, and same size. The objective function is generically defined as
%   the sum of the squared differences between the model and the tensor at
%   each element. The following parameter-value pairs define the objective
%   function and constraint details...   
%
%   o 'unique'    - In an M-way symmetric tensor, an element may appear up
%                   to M! times. This parameter controls whether or not to
%                   give each unique index equal weight. Otherwise, each  
%                   unique index is weighted by the number of appearances
%                   in the symmetric tensor. Default: True.
%   o 'fast'      - Use fast version of code if 'unique' is false. Not
%                   compatible with weights. Default: True.
%   o 'l2weight'  - Weight for the penalty term that is defined by
%                   sum_k (norm(X(:,k))^2 - 1)^2. Encourage column norms in
%                   X to be 1. Default: 0.  
%   o 'l1weight'  - Weight for the penalty term that is defined by 
%                   sum(LAMBDA). Encourages sparsity in LAMBDA. 
%                   Default: 0.
%   o 'l1param'   - Alpha-term in L1 approximation for penalizing
%                   encouraging sparsity in LAMBDA. Default: 10.
%   o 'nonneg'    - Specify nonnegativity constraints: X >= 0, LAMBDA >= 0.
%                   Default: False.
%   o 'nolambda'  - Remove LAMBDA from the optimization. Changes number of
%                   optimizations varaibles to Q = P*N. Default: False.
%   o 'weights'   - Specity weight tensor for the optimization. Must also
%                   be symmetric. Default: [].
%
%   Example:
%   A = symmetrize(tenrand(2,2,2)); P = 3; Model = symktensor(P,A);
%   data = fg_setup(MODEL,A);
%
%   See also SYMKTENSOR, CP_SYM, SYMKTENSOR/FG.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


   
% --- Process inputs ---
params = inputParser;
params.addParameter('unique', true, @islogical);
params.addParameter('fast', true, @islogical);
params.addParameter('l2weight', 0, @(x) isscalar(x) && x >= 0);
params.addParameter('l1weight', 0, @(x) isscalar(x) && x >= 0);
params.addParameter('l1param', 10, @(x) isscalar(x) && x >= 0);
params.addParameter('nonneg', false);
params.addParameter('nolambda', false, @islogical);
params.addParameter('weights',[]);
params.parse(varargin{:});
fgopts = params.Results;

% Check tensor
if ~(isa(A,'tensor')  && issymmetric(A)) && ~isa(A,'symtensor')
    error('A must be a symmetric tensor');
end

% Check weights
if ~isempty(fgopts.weights)
    if ~issymmetric(fgopts.weights)
        error('Weights must be symmetric');
    end
    if ~isequal(size(A),size(fgopts.weights))
        error('Weight tensor must be the same size as A');
    end
end

% Force fast option to be false if unique is true
if fgopts.unique
    fgopts.fast = false;
end

% Check fast option
if fgopts.fast
    if ~isempty(fgopts.weights)
        error('Cannot specify weights when ''fast'' is true');
    end
end

% Extract sizes
M = ndims(A);
N = size(A,1);
P = ncomponents(Model);

% Check lambda options
if fgopts.nolambda
    if (mod(M,2) == 0) && ~fgopts.nonneg
        error('Cannot have nolamdba=true and nonneg=false');
    end
    if fgopts.l2weight
        error('Cannot have nolambda=true and l2weight>0');
    end
    if fgopts.l1weight
        error('Cannot have nolambda=true and l1weight>0');
    end
end

if fgopts.nolambda
    R = N*P;
else
    R = (N+1)*P; % Number of free variables in model
end

% --- Bounds from constraints ---

% Variable bounds
lb = -Inf * ones(R,1);
ub = Inf * ones(R,1);

% Nonnegative
if (fgopts.nonneg)
    lb = zeros(R,1);
end

% --- Compute index sets and weights ---

if fgopts.fast
    
    % Convert to dense tensor
    A = full(A);
    
    % Compute squared norm
    normAsqr = norm(A)^2;
    
    % --- Assemble results & return ---
    data = var2struct(fgopts,M,N,P,A,normAsqr,lb,ub);

else
    A = symtensor(A);
    [I,C,W,Q] = indices(A);
    if fgopts.unique
        W = ones(size(W));
    end
       
    % --- Extract A values ---
    avals = A.val;
    
    % --- Incorporate weights if specified ---
    if ~isempty(fgopts.weights)
        if isa(fgopts.weights,'symtensor')
            W = W .* fgopts.weights.val;
        else
            W = W .* fgopts.weights(I);
        end
    end
    
    % --- Indices for gradient computation: onesidx and zeroidx ---
    % IDX(q,n) = mode of index that equal n in row q; 0 if DNE
    idx = zeros(Q,N);
    for n = 1:N
        [II,JJ] = find(I == n);
        idx(:,n) = accumarray(II,JJ,[Q 1],@min,0);
    end
    tf = idx > 0;
    
    len = zeros(N,1);
    idx1 = zeros(Q,N);
    idx2 = zeros(Q,N);
    for n = 1:N
        len(n) = sum(tf(:,n) > 0);
        idx1(1:len(n),n) = find(tf(:,n));
        idx1(len(n)+1:end,n) = find(~tf(:,n));
        idx2(1:len(n),n) = idx(tf(:,n),n);
    end
    
    onesidx = cell(n,1);
    zerosidx = cell(n,1);
    for n = 1:N
        onesidx{n} = tt_sub2ind([Q M],[idx1(1:len(n),n),idx2(1:len(n),n)]);
        zerosidx{n} = idx1(len(n)+1:end,n);
    end
    
    % --- Assemble results & return ---
    data = var2struct(fgopts,M,N,P,Q,R,avals,I,C,W,onesidx,zerosidx,lb,ub);
end

function s = var2struct(varargin)
%VAR2STRUCT
%http://stackoverflow.com/questions/3470654/how-can-i-move-variables-into-and-out-of-a-structure-akin-to-load-and-save-in-ma
names = arrayfun(@inputname,1:nargin,'UniformOutput',false);
s = cell2struct(varargin,names,2);

