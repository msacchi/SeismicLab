function t = uminus(t)
%UMINUS Unary minus for symktensor. 
%
%   See also SYMKTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



t.lambda = -t.lambda;
