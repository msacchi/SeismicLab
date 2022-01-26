function Z = tenfun(fun,varargin)
%TENFUN Apply a function to each element in a symmetric tensor.
%
%   TENFUN(F,X,...) applies the function specified by the function handle F
%   to the given arguments.  All arguments must be symtensors or scalars.
%   The functions are applied to the value arrays. If there are more than
%   two arguments, then the functions are applied elementwise.
%
%   Examples
%   Z = TENFUN(@(x)(x+1),X) %<-- increase every element by one
%   Z = TENFUN(@eq,X,1) %<-- logical comparison of X with scalar
%   Z = TENFUN(@plus,X,Y) %<-- adds the two symmetric tensors X and Y
%   Z = TENFUN(@max,X,Y,Z) %<-- max over all elements in X,Y,Z
%
%   See also SYMTENSOR, TENSOR/TENFUN.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if nargin < 2
    error('TTB:Symtensor:FunFail','Not enough input arguments.');
end

if ~isa(fun, 'function_handle')
    error('TTB:Symtensor:FunFail','First argument must be a function handle.');
end

%% Unary
if nargin == 2
    X = varargin{1}; % Must be a symtensor!
    Zvals = fun(X.val);
    Z = symtensor(Zvals, ndims(X), size(X,1));
    return;
end

%% Otherwise, sort out scalars and symtensors
tfscalar = cellfun(@isscalar, varargin);
tfsymtensor = cellfun(@(x) isa(x,'symtensor'), varargin);
sz = cellfun(@size, varargin(tfsymtensor), 'UniformOutput', false);

if ~all(tfscalar | tfsymtensor) 
    error('TTB:BadInputs','All inputs must be either symtensors or scalars');
end

if length(sz) > 1 && ~isequal(sz{:})
    error('TTB:BadInputs','All tensor inputs must be the same size');
end

m = length(sz{1});
n = sz{1}(1);

%% Binary Function
if nargin == 3
    X = cell(2,1);
    for j = 1:2
        if tfscalar(j)
            X{j} = varargin{j};
        else
            X{j} = varargin{j}.val;
        end
    end
    Zvals = fun(X{1},X{2});
    Z = symtensor(Zvals, m, n);
    return;
end

%% More than two inputs --- handle elementwise, if possible
p = nchoosek(m+n-1,m);
X = zeros(p,nargin-1);
for j = 1:nargin-1
    if tfscalar(j)
        X(:,j) = varargin{j};
    else
        X(:,j) = varargin{j}.val;
    end
end
Zvals = zeros(p,1);
for i = 1:p
    Zvals(i) = fun(X(i,:));
end
Z = symtensor(Zvals, m, n);
return;


