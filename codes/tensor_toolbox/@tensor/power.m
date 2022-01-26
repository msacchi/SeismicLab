function Z = power(X,Y)
%POWER Elementwise power (.^) operator for a tensor.
%
%   See also TENSOR.
% 
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



Z = tenfun(@power,X,Y);
