function Z = power(X,Y)
%POWER Elementwise power (.^) operator for a symmetric tensor.
%
%   See also SYMTENSOR.
% 
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


Z = tenfun(@power,X,Y);
