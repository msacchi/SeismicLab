function X = tendiag(v,sz)
%TENDIAG Creates a tensor with v on the diagonal.
%
%   TENDIAG(V) creates a tensor with N dimensions, each of size N, where N
%   is the number of elements of V. The elements of V are placed on the
%   superdiagonal.
%
%   TENDIAG(V,SZ) is the same as above but creates a tensor of size SZ. If
%   SZ is not big enough, the tensor will be enlarged to accommodate the
%   elements of V on the superdiagonal.
%
%   Examples
%   X = tendiag([0.1 0.22 0.333]) %<-- creates a 3x3x3 tensor
%
%   See also TENSOR, SPTENDIAG.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


% Make sure v is a column vector
v = reshape(v,[numel(v) 1]);

N = numel(v);
if ~exist('sz','var')
    sz = repmat(N,1,N);
end

X = tenzeros(sz);
subs = repmat((1:N)', 1, length(sz));
X(subs) = v;
