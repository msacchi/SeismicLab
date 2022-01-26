function sz = size(t,idx)
%SIZE Size of symktensor.
%  
%   D = SIZE(T) returns the size of the tensor. 
%
%   I = SIZE(T,DIM) returns the size of the dimension specified by
%   the scalar DIM.
%
%   See also SYMKTENSOR, SYMKTENSOR/NDIMS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isempty(t.lambda)
    sz = [];
end

if exist('idx','var')
    sz = size(t.u, 1);
else
    sz = size(t.u, 1) * ones(1,t.m);  
end
