function n = ndims(t)
%NDIMS Number of dimensions of a sparse tensor.
%
%   NDIMS(T) returns the number of dimensions of sparse tensor T.  
%
%   Examples:
%   T = sptenrand([3 2 2],5); 
%   ndims(T) %<-- should return 3
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



n = size(t.size,2);
