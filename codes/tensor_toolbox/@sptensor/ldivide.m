function C = ldivide(A,B)
%LDIVIDE Array right division for sparse tensors.
%
%   LDIVIDE(A,B) is called for the syntax 'A .\ B' when A or B is a sparse
%   tensor. A and B must have the same size, unless one is a scalar. 
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



C = rdivide(B,A);
