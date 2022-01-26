function n = norm(T)
%NORM Frobenius norm of a sumtensor.
%
%   It is not possible to efficiently compute the NORM of sumtensor. We
%   therefore just return zero and print a warning. This function is
%   included for compatibility with certain routines that expect to compute
%   the norm but don't *really* need it.
%   
%   NORM(X) returns 0.
%
%   See also SUMTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



warning('The NORM function is not supported by SUMTENSOR. Returning zero.');
n = 0;
