function siz = size(a,idx)
%SIZE Return size of sptenmat.
%  
%   D = SIZE(T) returns the size of the tensor.  
%
%   I = size(T,DIM) returns the sizes of the dimensions specified by DIM.
%
%   See also SPTENMAT.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isempty(a.tsize)
    siz = [];
    return;
end

m = prod(a.tsize(a.rdims));
n = prod(a.tsize(a.cdims));
siz = [m n];

if  exist('idx','var')
    siz = siz(idx);
end
