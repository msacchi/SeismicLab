function n = ndims(t)
%NDIMS Return the number of dimensions for a sumtensor.
%
%   NDIMS(T) returns the number of dimensions of tensor T.
%
%   See also SUMTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


n = ndims(t.part{1});
