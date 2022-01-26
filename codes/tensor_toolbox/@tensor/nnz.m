function n = nnz(x)
%NNZ Number of nonzeros for tensors. 
%
%   See also TENSOR, SPTENSOR/NNZ.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



n = nnz(x.data);
