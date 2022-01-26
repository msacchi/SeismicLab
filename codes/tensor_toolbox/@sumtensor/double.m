function A = double(T)
%DOUBLE Convert sumtensor to double array.
%
%   A = double(T) converts T to a standard multidimensional array. 
%
%   See also SUMTENSOR, SUMTENSOR/FULL.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


A = double(full(T));
