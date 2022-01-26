function X = tenzeros(varargin)
%TENZEROS Create zeros tensor.
%
%   X = TENZEROS(SZ) forms a tensor of size SZ with all zeros.
%
%   TENZEROS(SZ) is equivalent to TENSOR(ZEROS(SZ(1),SZ(2),...),SZ).
%
%   See also TENSOR, ZEROS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if nargin == 1
    sz = varargin{1};
else
    sz = cell2mat(varargin);
end

if isempty(sz)
    X = tensor;
    return;
end

if nargin == 2
    order = sz;
    dim = varargin{1};
    sz = dim * ones(1,order);
end

data = zeros([sz 1 1]);
X = tensor(data,sz);

