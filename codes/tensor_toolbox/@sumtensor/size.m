function m = size(t,idx)
%SIZE Size of a sumtensor.
%  
%   D = SIZE(T) returns the size of the tensor. 
%
%   I = size(T,DIM) returns the size of the dimension specified by
%   the scalar DIM.
%
%   See also SUMTENSOR, SUMTENSOR/NDIMS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if isempty(t.part)
    m = [];
elseif exist('idx','var')
    m = size(t.part{1},idx);
else
    m = size(t.part{1});
end
    
