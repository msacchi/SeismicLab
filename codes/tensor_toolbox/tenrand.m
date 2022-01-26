function X = tenrand(varargin)
%TENRAND Uniformly distributed pseudo-random tensor.
%
%   X = TENRAND(SZ) forms a tensor of size SZ with pseudo-random
%   values drawn from a uniform distribution on the unit interval.
%
%   TENRAND(SZ) is equivalent to TENSOR(RAND(SZ(1),SZ(2),...),SZ).
%
%   See also TENSOR, SPTENRAND, RAND.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if nargin == 1
    sz = varargin{1};
else
    sz = cell2mat(varargin);
end

if isempty(sz)
    X = tensor;
else
    data = rand([sz 1 1]);
    X = tensor(data,sz);
end
