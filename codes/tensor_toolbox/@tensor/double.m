function A = double(X)
%DOUBLE Convert tensor to double array.
%
%   A = double(X) converts X to a standard multidimensional array.
%
%   See also TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



A = double(X.data);
