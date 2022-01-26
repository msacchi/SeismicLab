%SYMKTENSOR Class for storing symmetric Kruskal tensor (decomposed).
%
%SYMKTENSOR Methods:
%   arrange     - Arranges the rank-1 components of a symktensor.
%   disp        - Command window display for a symktensor.
%   display     - Command window display for a symktensor.
%   double      - Convert a symktensor to a double array.
%   end         - Last index of indexing expression for symktensor.
%   entry       - Extract a single entry from a symktensor.
%   f_implicit  - Function F(M)=||X-M||^2 for two symktensors.
%   fg          - Master objective function for optimization of symmetric Kruskal model.
%   fg_explicit - Function and gradient of F(M)=||X-M||^2 for symktensor.
%   fg_implicit - Function and gradient of F(M)=||X-M||^2 for two symktensors.
%   fg_setup    - Setup for optimization of symmetric Kruskal model.
%   full        - Convert a symktensor to a symtensor.
%   g_implicit  - Gradient of F(M)=||X-M||^2 for two symktensors.
%   isequal     - True if each component of two symktensors is numerically equal.
%   isscalar    - False for symktensors.
%   issymmetric - Rhetorical function for a symktensor.
%   mtimes      - Implement A*B (scalar multiply) for symktensor.
%   ncomponents - Number of components for a symktensor.
%   ndims       - Number of modes for a symktensor.
%   norm        - Frobenius norm of a symktensor.
%   normalize   - Normalizes the columns of the factor matrix.
%   permute     - Permute dimensions of a symktensor.
%   randextract - Create new symktensor with random subset of components.
%   score       - Checks if two symktensors match except for permutation.
%   size        - Size of symktensor.
%   subsasgn    - Subscripted assignment for symktensor.
%   subsref     - Subscripted reference for a symktensor.
%   symktensor  - Tensor stored as a symmetric Kruskal operator (decomposed).
%   tovec       - Convert symktensor to vector representation.
%   uminus      - Unary minus for symktensor. 
%   update      - Convert vector to symktensor.
%   uplus       - Unary plus for a symktensor. 
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','symktensor_doc.html')))">Documentation page for symmetric Kruskal tensor class</a>
%
%   See also TENSOR_TOOLBOX, SYMTENSOR
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

function t = symktensor(varargin)
%SYMKTENSOR Tensor stored as a symmetric Kruskal operator (decomposed).
%
%   A symmetric Kruskal tensor is used to build a model of a M-way
%   N-dimensional symmetric tensor. We have to specify the number of
%   components, denoted by P. Each component comprises a scalar weight and
%   the M-way symmetric outer product of a vector. We store the P weights
%   in the length-P vector LAMBDA and the vectors in the NxP matrix X. The
%   number of modes, M, is stored as a scalar.
%
%   S = SYMKTENSOR(LAMBDA,U,M) creates a symmetric M-way Kruskal tensor
%   from its constituent parts. Here lambda is a K-vector and each U is a  
%   matrix with K columns.
%
%   S = SYMKTENSOR(K) creates a symmetric Kruskal tensor by "symmetrizing"
%   a ktensor, i.e., averaging the constituent factor matrices and taking
%   care to get the signs aligned. 
%
%   S = SYMKTENSOR(S0) creates a symktensor by copying an existing
%   symktensor. 
%
%   *** Below are specialized constructore for use in optimization. There
%   are not recommended for general use. ***
%
%   S = SYMKTENSOR(P,A) creates a symmetric Krukal tensor with P components
%   that is sized to match the tensor A. The LAMBDA is set to be all ones,
%   and the X is initialized using to uniform random values in [0,1].
%
%   S = SYMKTENSOR(V,A) creates a symmetric Kruskal tensor from a
%   vectorized verion V. The second argument is a tensor of the size that
%   is being modeled.
%
%   S = SYMKTENSOR(V,A,NOLAMBDA) creates a symmetric Kruskal tensor from a
%   vectorized version V. The second argument is a tensor of the size that
%   is being modeled. The third argument indicates whether or not LAMBDA is
%   stored in V. If NOLAMBDA=TRUE, then only X is stored (i.e., LAMBDA
%   values are set to 1).  
%
%   S = SYMKTENSOR(V,M,P) creates an M-way symmetric Kruskal
%   tensor with P components from a vectorized version V.  
%
%   S = SYMKTENSOR(V,M,P,NOLAMBDA) uses the NOLAMBDA parameter as described
%   above. 
%
%   Examples
%   LAMBDA = randn(2, 1); %<-- Lambda vector (should be a column)
%   U = randn(4, 2); %<-- Factor matrix
%   S = symktensor(LAMDBA, U, 3) %<--Declare symktensor 
%
%   See also KTENSOR, TOVEC
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



% EMPTY CONSTRUCTOR
if nargin == 0
    t.lambda = [];
    t.u = [];
    t.m = 0;
    t = class(t,'symktensor');
    return;
end

% Copy CONSTRUCTOR
if (nargin == 1) && isa(varargin{1}, 'symktensor')
    t.lambda = varargin{1}.lambda;
    t.u = varargin{1}.u;
    t.m = varargin{1}.m;
    t = class(t, 'symktensor');
    return;
end

% Symmetrize a KTENSOR
if (nargin == 1) && isa(varargin{1}, 'ktensor')
    K = varargin{1};
    if ~issymmetric(K)
         K = symmetrize(K);
    end
    t.lambda = K.lambda;
    t.u = K.U{1};
    t.m = length(K.U);
    t = class(t, 'symktensor');
    return;
end

% Create a random SYMKTENSOR from the number of components and a tensor.
if (nargin == 2) && isscalar(varargin{1})
    P = varargin{1};
    A = varargin{2};
    if  ~isa(A,'tensor') && ~isa(A,'symtensor')
        error('Second argument should be a tensor when the first argument is a scalar');
    end
    M = ndims(A);
    N = size(A,1);
    t.lambda = ones(P,1); %Column vector
    t.u = rand(N,P);   %Robert development note: Document this better
    t.m = M;
    t = class(t,'symktensor');
    return;
end

% Convert from vector representation, given tensor to determine shape
if nargin == 2
   v = varargin{1};
   A = varargin{2}; 
   if ~isvector(v) || ~(isa(varargin{2},'tensor') || isa(varargin{2},'symtensor'))
       error('Wrong arguments');
   end
   n = size(A,1);
   p = length(v) / (n+1);
   if round(p) ~= p
       error('Length of input vector is not evenly divisible by (n+1)');
   end
   t.lambda = v(1:p); 
   u = v(p+1:end);
   t.u = reshape(u,n,p);
   t.m = ndims(A);
   t = class(t, 'symktensor');
   return;
end


% Convert from vector representation, given tensor to determine shape
if (nargin == 3) && (isa(varargin{2},'tensor') || isa(varargin{2},'symtensor'))
    
    tf = varargin{3}; %No lambda option
    
    if ~tf
        t = symktensor(varargin{1},varargin{2});
        return;
    end
    
   v = varargin{1};
   A = varargin{2};
   if ~isvector(v)
       error('Wrong arguments');
   end
   n = size(A,1);
   p = length(v) / n;
   if round(p) ~= p
       error('Length of input vector is not evenly divisible by n');
   end
   t.lambda = ones(p,1); %Column vector
   u = v;
   t.u = reshape(u,n,p);
   t.m = ndims(A);
   t = class(t, 'symktensor');
   return;
end

% Convert from vector representation, given sizes
if (nargin == 3) && isscalar(varargin{2})
    v = varargin{1};
    m = varargin{2};
    p = varargin{3};
    t.lambda = v(1:p); 
    u = v(p+1:end);
    t.u = reshape(u,[],p);
    t.m = m;
    t = class(t, 'symktensor');
    return;
end

% Convert from vector representation, given sizes and nolambda!
if (nargin == 4) 
    tf = varargin{4};
    
    if ~tf
        t = symktensor(varargin{1},varargin{2},varargin{3});
        return;
    end
    
    v = varargin{1};
    m = varargin{2};
    p = varargin{3};
    t.lambda = ones(p,1); %Column vector
    u = v;
    t.u = reshape(u,[],p);
    t.m = m;
    t = class(t, 'symktensor');
    return;
end


if nargin ~= 3
    error('Check arguments to create symktensor');
end

t.lambda = varargin{1};
t.u = varargin{2};
t.m = varargin{3};

if ~isa(t.lambda,'numeric') || ~iscolumn(t.lambda)
    error('LAMBDA must be a column vector.'); 
end
    
% Check that each Um is indeed a matrix
if ~ismatrix(t.u) 
	error(['Matrix U is not a matrix!']);
end

% Size error checking			     
k = length(t.lambda); 
if  size(t.u,2) ~= k
       error(['Matrix U does not have ' int2str(k) ' columns.']);
end

t = class(t, 'symktensor');
return;
