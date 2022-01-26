function t = uminus(t)
%UMINUS Unary minus (-) for sptensor.
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



t.vals = -t.vals;
