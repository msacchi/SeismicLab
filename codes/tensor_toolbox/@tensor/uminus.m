function t = uminus(t)
%UMINUS Unary minus (-) for tensors.
%
%   See also TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



t.data = -t.data;
