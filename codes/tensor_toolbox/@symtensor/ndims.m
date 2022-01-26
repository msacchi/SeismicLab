function n = ndims(t)
%NDIMS Number of dimensions for a symtensor.
%
%   NDIMS(S) returns the number of dimensions of a symtensor S.
%
%   Examples
%   X = symtensor(rand([4,4,4])); ndims(X) %<-- Returns 3
%
%   See also SYMTENSOR, TENSOR/NDIMS
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



n = t.m;
