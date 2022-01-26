function nrm = norm(T)
%NORM Frobenius norm of a sparse tensor.
%
%   NORM(T) returns the Frobenius norm of a sparse tensor.
%
%   See also SPTENSOR, NORM.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



nrm = norm(T.vals);

return;
