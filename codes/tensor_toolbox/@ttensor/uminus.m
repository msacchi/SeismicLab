function t = uminus(t)
%UMINUS Unary minus for ttensor.
%
%   See also TTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


t.core = -t.core;
