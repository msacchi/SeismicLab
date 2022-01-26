function A = double(T)
%DOUBLE Convert ttensor to double array.
%
%   A = double(T) converts T to a standard multidimensional array.
%
%   See also TTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


A = double(full(T));
