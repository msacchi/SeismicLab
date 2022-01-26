function Z = tenfun(fun,varargin)
%TENFUN Apply a function to each element in a tensor.
%
%   TENFUN(F,X,...) applies the function specified by the function
%   handle F to the given arguments.  Either both arguments
%   must be tensors, or one is a tensor and the other is a scalar/MDA.
%
%   Examples
%   Z = tenfun(@(x)(x+1),X) %<-- increase every element by one
%   Z = tenfun(@eq,X,1) %<-- logical comparison of X with scalar
%   Z = tenfun(@plus,X,Y) %<-- adds the two tensors X and Y.
%   Z = tenfun(@max,X,Y,Z) %<-- elementwise max over all elements in X,Y,Z
%
%   See also TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if nargin < 2
    %  error('TENFUN requires at least two input arguments')
    error('Not enough input arguments.');
end

if ~isa(fun, 'function_handle')
    error('First argument must be a function handle.');
end

%% Convert inputs to tensors if they aren't already 
for i = 1:nargin-1
    if isscalar(varargin{i}) || isa(varargin{i},'tensor')
        continue;
    elseif isnumeric(varargin{i})  
        varargin{i} = tensor(varargin{i});
    elseif ismember(class(varargin{i}), {'symtensor','sptensor','ktensor','ttensor'})
        varargin{i} = full(varargin{i});
    else
        error('Invalid input');        
    end        
end
%% It's ok if there are two inputs and one is a scalar; otherwise, all inputs have to be the same size
if (nargin == 3) && isscalar(varargin{1}) && isa(varargin{2},'tensor')
    sz = size(varargin{2});
elseif (nargin == 3) && isscalar(varargin{2}) && isa(varargin{1},'tensor')
    sz = size(varargin{1});
else
    for i = 1:(nargin-1)
        if isscalar(varargin{i})
            error('Argument %d is a scalar, but expected a tensor', i+1);
        elseif i == 1
            sz = size(varargin{i});
        elseif ~isequal(sz,size(varargin{i}))
            error('Tensor %d is not the same size as the first tensor input', i);
        end
    end
end

%% Number of inputs for function handle
nfunin = nargin(fun);

%% Case I: Binary function
if (nargin == 3) && (nfunin == 2)   
    X = varargin{1};
    Y = varargin{2};    
    if ~isscalar(X)
        X = X.data;
    end
    if ~isscalar(Y)
        Y = Y.data;
    end
    data = fun(X,Y);
    Z = tensor(data, sz);
    return;
end


%% Case II: Expects input to be matrix and applies operation on each column
if nargin == 2
    X = varargin{1}.data;
    X = reshape(X,1,[]);
else
    X = zeros(nargin-1,prod(sz));
    for i = 1:nargin-1
        X(i,:) = varargin{i}.data(:);
    end
end
data = fun(X);
Z = tensor(data,sz);
return;



