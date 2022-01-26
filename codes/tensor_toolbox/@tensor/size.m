function m = size(t,idx)
%SIZE Tensor dimensions.
%  
%   D = SIZE(T) returns the sizes of each dimension of tensor X in a
%   vector D with ndims(X) elements.
%
%   I = size(T,DIM) returns the size of the dimension specified by
%   the scalar DIM.
%
%   Examples
%   A = rand(3,4,2,1); T = tensor(A,[3 4 2 1]);
%   size(A) %<-- returns a length-3 vector
%   size(T) %<-- returns a length-4 vector
%   size(A,2) %<-- returns 4
%   size(T,2) %<-- same
%   size(A,5) %<-- returns 1
%   size(T,5) %<-- ERROR!
%
%   See also TENSOR, TENSOR/NDIMS, SIZE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>




if exist('idx','var')
    m = t.size(idx);
else
    m = t.size;
end
