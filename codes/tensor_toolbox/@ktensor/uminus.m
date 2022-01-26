function t = uminus(t)
%UMINUS Unary minus for ktensor. 
%
%   See also KTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



t.lambda = -t.lambda;
