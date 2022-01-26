function n = nnz(a)
%NNZ Return number of nonzeros in a sptenmat.
%
%   nnz(A) returns the number of nonzeros in A.
%
%   See also SPTENMAT.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



n = length(a.vals);
