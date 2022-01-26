function X = tenones(varargin)
%TENONES Ones tensor.
%
%   X = TENONES(SZ) forms a tensor of size SZ with all ones.
%
%   TENONES(SZ) is equivalent to TENSOR(ONES(SZ(1),SZ(2),...),SZ).
%
%   See also TENSOR, ONES.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if nargin == 1
    sz = varargin{1};
else
    sz = cell2mat(varargin);
end

if isempty(sz)
    X = tensor();
    return;
end

data = ones([sz 1 1]);
X = tensor(data,sz);
