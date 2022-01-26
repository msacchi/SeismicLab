function t = exp(t)
%EXP Exponential for tensors.
%
%   EXP(X) is the exponential of the elements of X.
%
%   See also TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


t.data = exp(t.data);
