function n = ndims(t)
%NDIMS Number of dimensions for a ktensor.
%
%   NDIMS(T) returns the number of dimensions of tensor T.
%
%   See also KTENSOR
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



n = numel(t.u);
