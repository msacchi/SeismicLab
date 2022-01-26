function nrm = norm(T)
%NORM Frobenius norm of a sptenmat.
%
%   NORM(T) returns the Frobenius norm of a matricized sparse tensor.
%
%   See also SPTENMAT, NORM.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



nrm = norm(T.vals);

return;
