function n = ndims(t)
%NDIMS Return the number of dimensions of a tensor.
%
%   NDIMS(X) returns the number of dimensions of tensor X.
%
%   Examples
%   A = rand(4,3,1); ndims(A) %<-- Returns 2
%   X = tensor(A); ndims(X) %<-- Returns 2
%   X = tensor(A,[4 3 1]); ndims(X) %<-- Returns 3
%
%   See also TENSOR
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



n = numel(t.size);
