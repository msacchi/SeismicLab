function sz = tsize(a,idx)
%TSIZE Tensor size of sptenmat.
%
%   D = TSIZE(X) returns the size of the tensor being stored as a
%   matrix. 
% 
%   M = TSIZE(X,DIM) returns the length of the dimension(s) specified
%   by DIM.  For example, SIZE(X,1) returns the size of the first
%   dimension of the tensor.
%
%   See also SPTENMAT, SPTENMAT/SIZE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isempty(a.tsize)
    sz = [];
    return;
end

if exist('idx', 'var')
    sz = a.tsize(idx);
else
    sz = a.tsize;
end

return;
