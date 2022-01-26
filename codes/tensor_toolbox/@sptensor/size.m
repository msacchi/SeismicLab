function m = size(t,idx)
%SIZE Sparse tensor dimensions.
%  
%   D = SIZE(T) returns the size of the tensor.  
%
%   I = size(T,DIM) returns the sizes of the dimensions specified by DIM,
%   which is either a scalar or a vector of dimensions.
%
%   See also SPTENSOR, SPTENSOR/NDIMS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if exist('idx','var')
    m = t.size(idx);
else
    m = t.size;
end
