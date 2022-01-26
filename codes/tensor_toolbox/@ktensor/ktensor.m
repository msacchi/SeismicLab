%KTENSOR Class for Kruskal tensors (decomposed).
%
%KTENSOR Methods:
%   arrange      - Arranges the rank-1 components of a ktensor.
%   datadisp     - Special display of a ktensor.
%   disp         - Command window display for a ktensor.
%   display      - Command window display for a ktensor.
%   double       - Convert a ktensor to a double array.
%   end          - Last index of indexing expression for ktensor.
%   extract      - Creates a new ktensor with only the specified components.
%   fixsigns     - Fix sign ambiguity of a ktensor.
%   full         - Convert a ktensor to a (dense) tensor.
%   innerprod    - Efficient inner product with a ktensor.
%   isequal      - True if each datum of two ktensor's are numerically equal.
%   isscalar     - False for ktensors.
%   issymmetric  - Verify that a ktensor X is symmetric in all modes.
%   ktensor      - Tensor stored as a Kruskal operator (decomposed).
%   mask         - Extract values as specified by a mask tensor.
%   minus        - Binary subtraction for ktensor.  
%   mtimes       - Implement A*B (scalar multiply) for ktensor.
%   mttkrp       - Matricized tensor times Khatri-Rao product for ktensor.
%   ncomponents  - Number of components for a ktensor.
%   ndims        - Number of dimensions for a ktensor.
%   norm         - Frobenius norm of a ktensor.
%   normalize    - Normalizes the columns of the factor matrices.
%   nvecs        - Compute the leading mode-n vectors for a ktensor.
%   permute      - Permute dimensions of a ktensor.
%   plus         - Binary addition for ktensor.
%   redistribute - Distribute lambda values to a specified mode. 
%   score        - Checks if two ktensors match except for permutation.
%   size         - Size of ktensor.
%   subsasgn     - Subscripted assignment for ktensor.
%   subsref      - Subscripted reference for a ktensor.
%   symmetrize   - Symmetrize a ktensor X in all modes.
%   times        - Element-wise multiplication for ktensor.
%   tocell       - Convert X to a cell array.
%   tovec        - Convert Ktensor to vector.
%   ttm          - Tensor times matrix for ktensor.
%   ttv          - Tensor times vector for ktensor.
%   uminus       - Unary minus for ktensor. 
%   update       - Update one or more modes of the ktensor with new data.
%   uplus        - Unary plus for a ktensor. 
%   viz          - Visualize a ktensor.
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','ktensor_doc.html')))">Documentation page for Kruskal tensor class</a>
%
%   See also TENSOR_TOOLBOX
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

function t = ktensor(varargin)
%KTENSOR Tensor stored as a Kruskal operator (decomposed).
%
%   K = KTENSOR(lambda,U1,U2,...,UM) creates a Kruskal tensor from its
%   constituent parts. Here lambda is a k-vector and each Um is a
%   matrix with k columns.
%
%   K = KTENSOR(lambda, U) is the same as above except that U is a
%   cell array containing matrix Um in cell m.
%
%   K = KTENSOR(U) assumes U is a cell array containing matrix Um in
%   cell m and assigns the weight of each factor to be one.
%
%   K = KTENSOR(T) creates a ktensor by copying an existing ktensor.
%
%   K = KTENSOR(S) creates a ktensor from a symktensor.
%
%   K = KTENSOR(FH, SZ, NC) creates a ktensor using function FH to create
%   the factor matrices. Here SZ is the size of the final ktensor and NC is
%   the number of components. The function specified by FH should take two
%   size arguments and create a matrix of that size.
%
%   Examples
%   K = ktensor([3; 2], ones(4,2), ones(5,2), ones(3,2)) %<- Constructor
%   K = ktensor(@rand, [4 5 3], 2) %<- Create a random tensor
%
%   See also KTENSOR, CP_ALS, CP_OPT, CP_WOPT, CP_APR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



% EMPTY CONSTRUCTOR
if nargin == 0
    t.lambda = [];
    t.u = {};
    t = class(t,'ktensor');
    return;
end

% Copy CONSTRUCTOR
if (nargin == 1) && isa(varargin{1}, 'ktensor')
    t.lambda = varargin{1}.lambda;
    t.u = varargin{1}.u;
    t = class(t, 'ktensor');
    return;
end


% CONSTRUCTOR from SYMKTENSOR
% TODO: Need to check that this works! 
if (nargin == 1) && isa(varargin{1}, 'symktensor')
    t.lambda = varargin{1}.lambda;
    [t.u{1:varargin{1}.m,1}] = deal(varargin{1}.u);
    t = class(t, 'ktensor');
    return;
end

% CONSTRUCTOR FROM FUNCTION HANDLE
if (nargin == 3) && isa(varargin{1},'function_handle')
    fh = varargin{1};
    sz = varargin{2};
    nc = varargin{3};
    nd = length(sz);
    lambda = ones(nc,1);
    U = cell(nd,1);
    for i = 1:length(sz)
        U{i} = feval(fh,sz(i),nc);
    end
    t.lambda = lambda;
    t.u = U;
    t = class(t,'ktensor');
    return;
end

% CONSTRUCTOR from VECTOR, SIZE, NDIMS, LAMBDAFLAG
if (nargin == 4) && isvector(varargin{1}) && islogical(varargin{4})
    x = varargin{1};
    sz = varargin{2};
    nd = varargin{3};
    lambdaflag = varargin{4};        
    
    if isrow(x)
        x = x';
    end
    
    if lambdaflag 
        nc = length(x) / (sum(sz) + 1);
    else
        nc = length(x) / sum(sz);
    end
    if round(nc) ~= nc
        error('Vector is not the right length');
    end
    
    if lambdaflag
        t.lambda = x(1:nc);
        shift = nc;
    else
        t.lambda = ones(nc,1);
        shift = 0;
    end

    t.u = cell(nd,1);
    for n = 1:nd
        mstart = nc*sum(sz(1:n-1))+shift+1;
        mend = nc*sum(sz(1:n))+shift;
        t.u{n} = reshape(x(mstart:mend),[],nc);
    end
    t = class(t, 'ktensor');
    return;
end

% CONSTRUCTOR FROM CELL ARRAY
if (nargin == 1) && isa(varargin{1},'cell')

    u = varargin{1};
    nc = size(u{1},2);
    if ~all(cellfun(@(x) ismatrix(x) && size(x,2) == nc, u))
        error('Invalid factor matrix')
    end
    t.lambda = ones(nc,1);
    t.u = u;
    if ~isvector(t.u)
        error('U must be a vector');
    end
    if isrow(t.u)
        t.u = t.u';
    end
    t = class(t, 'ktensor');
    return;
    
end

% CONSTRUCTOR FOR LAMBDA and LIST OF MATRICES
if (nargin >= 2)

    t.lambda = varargin{1};
    if ~iscolumn(t.lambda)
        error('LAMBDA must be a column vector.');
    end
    
    if isa(varargin{2},'cell')
        t.u = varargin{2};
    else
        t.u = varargin(2:end);
    end
    
    if ~isvector(t.u)
        error('U must be a vector');
    end
    if isrow(t.u)
        t.u = t.u';
    end
    
    nc = length(t.lambda);
    if ~all(cellfun(@(x) ismatrix(x) && size(x,2) == nc, t.u))
        error('Invalid factor matrix')
    end

    t = class(t, 'ktensor');

    return;
end

error('Invalid ktensor constructor');
