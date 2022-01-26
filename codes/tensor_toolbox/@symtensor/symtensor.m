%SYMTENSOR Class for storing only unique entries of symmetric tensor.
%
%SYMTENSOR Methods:
%   and         - Logical AND (&) for symmetric tensors.
%   disp        - Command window display for a symtensor.
%   display     - Command window display of a symtensor.
%   eq          - Equal (==) for tensors.
%   full        - Convert symtensor to a tensor.
%   ge          - Greater than or equal (>=) for tensors.
%   gt          - Greater than (>) for symmetric tensors.
%   indices     - Compute unique indices of a symmetric tensor.
%   isequal     - for symmetric tensors.
%   isscalar    - False for symtensors.
%   issymmetric - Checks if tensor is symmetric (always true for symtensor).
%   ldivide     - Left array divide (.\) for symmetric tensors.
%   le          - Less than or equal (<=) for symmetric tensors.
%   lt          - Less than (<) for symmetric tensor.
%   minus       - Binary subtraction (-) for symmetric tensors.
%   mldivide    - Slash left division for symmetric tensors.
%   mrdivide    - Slash right division for symmetric tensors.
%   mtimes      - tensor-scalar multiplication.
%   ndims       - Number of dimensions for a symtensor.
%   ne          - Not equal (~=) for symmetric tensors.
%   not         - Logical NOT (~) for symmetric tensors.
%   or          - Logical OR (|) for symmetric tensors.
%   plus        - Binary addition (+) for symtensors.
%   power       - Elementwise power (.^) operator for a symmetric tensor.
%   rdivide     - Right array divide (./) for tensors.
%   size        - Dimensions of a symmetric tensor.
%   subsasgn    - Subassignment for symtensor.
%   subsref     - Subreference function for symtensor.
%   symtensor   - Symmetric tensor that stores only the unique values.
%   tenfun      - Apply a function to each element in a symmetric tensor.
%   times       - Array multiplication (.*) for symmetric tensors.
%   uminus      - Unary minus (-) for tensors.
%   uplus       - Unary plus (+) for symmetric tensors.
%   xor         - Logical EXCLUSIVE OR for symmetric tensors.
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','symtensor_doc.html')))">Documentation page for symmetric tensor class</a>
%
%   See also TENSOR_TOOLBOX, TENSOR/SYMMETRIZE, SYMKTENSOR, CP_SYM
%
%MATLAB Tensor Toolbox. Copyright 2017, Sandia Corporation.

function [t,I] = symtensor(varargin)
%SYMTENSOR Symmetric tensor that stores only the unique values.
%
%   S = SYMTENSOR(X) creates a symmetric tensor from a given tensor X. If X
%   is not symmetric, than it is symmetrized.
%
%   [S,I] = SYMTENSOR(X) when X is a tensor also returns I, a matrix of
%   indices of X which define the symmetric tensor. (This is the only case
%   where two arguments are returned.)
%
%   S = SYMTENSOR(S0) copies the symtensor S0.
%
%   S = SYMTENSOR(VALS,M,N) constructs a symmetric tensor with M modes,
%   dimension N, and values specified by VALS.
%
%   S = SYMTENSOR(FUN,M,N) constructs a symmetric tensor with M
%   modes, dimension N, and values generated with function FUN. The
%   function fun must return a matrix of values given arguments of row and
%   column size.
%
%   Examples
%   X = tenrand([3 3 3]);
%   S = symtensor(X); % Symmetrize X
%   S = symtensor(1:15,4,3); % Specify unique values, order 3, dim 3
%   S = symtensor(@rand,3,7); % Random, order 3, dim 7
%   S = symtensor(@ones,4,7); % Ones of order 4, dim 7
%
%   See also TENSOR/SYMMETRIZE, SYMKTENSOR, CP_SYM
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


% Set up --- the same for all
t = init_fields;
t = class(t,'symtensor');

% EMPTY CONSTRUCTOR
if nargin == 0
    error('TTB:BadInputs', 'Not enough input arguments.');
end

% COPY CONSTRUCTOR
if nargin == 1 && isa(varargin{1}, 'symtensor')
    t.val = varargin{1}.val;
    t.m = varargin{1}.m;
    t.n = varargin{1}.n;
    return;
end

% CREATE FROM TENSOR
if nargin == 1 && isa(varargin{1}, 'tensor')
    
    src = varargin{1};
    src = symmetrize(src);
    
    m = ndims(src);
    n = size(src,1);
    t.m = m;
    t.n = n;
    
    % Generate distinct indices
    sz = nchoosek(m+n-1,m);
    I = zeros(sz,m);
    for loc = 1:sz
        if loc == 1
            I(loc,:) = ones(1,m);
        else
            I(loc,:) = I(loc-1,:);
            j = m;
            while (I(loc,j) == n)
                j = j - 1;
            end
            I(loc,j:m) = I(loc,j)+1;
        end
    end
    
    % Query symmetric indices from tensor, then save into val array
    t.val = double(src(I));
    return;
end

% For constructing a symtensor from a function handle
if nargin==3 && isa(varargin{1},'function_handle')
    m = varargin{2};
    n = varargin{3};
    t.m = m;
    t.n = n;
    sz = nchoosek(m+n-1,m);
    try
        t.val = feval(varargin{1}, sz, 1);
        t.val = double(t.val); % Convert to double, if possible
    catch
        error('TTB:BadInputs','Bad generating function');
    end
    if ~isequal(size(t.val),[nchoosek(n+m-1,m),1])
        error('TTB:BadInputs','Bad generating function');
    end
    return;
end

% For constructing a symtensor from a value array
if nargin==3
    val = varargin{1};
    m = varargin{2};
    n = varargin{3};
    sz = nchoosek(m+n-1,m);
    if ~isvector(val) || numel(val) ~= sz
        error('TTB:BadInputs', 'Value array is the wrong size');
    end
    t.m = m;
    t.n = n;
    if iscolumn(val) % Needs to be a column array
        t.val = val;
    else
        t.val = val';
    end
    return;
end

error('TTB:BadInputs','Too many input arguments');
end

function t = init_fields()
t.val = [];
t.m = 0;
t.n = 0;
end


