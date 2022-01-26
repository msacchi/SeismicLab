function n = ndims(t)
%NDIMS Return the number of dimensions for a ttensor.
%
%   NDIMS(T) returns the number of dimensions of tensor T.
%
%   See also TTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


n = numel(t.u);
