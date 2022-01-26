function t = ones(t)
%ONES Replace nonzero elements of sparse tensor with ones.
%
%   S = ONES(T) generates a sparse tensor with the same sparsity
%   structure as T, but with ones in the nonzero position.
%
%   See also SPTENSOR, SPONES.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



t.vals = ones(size(t.vals));
