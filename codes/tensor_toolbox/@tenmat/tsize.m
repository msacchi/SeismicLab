function sz = tsize(a,idx)
%TSIZE Tensor size of tenmat.
%
%   D = TSIZE(X) returns the size of the tensor being stored as a
%   matrix. 
% 
%   M = TSIZE(X,DIM) returns the length of the dimension(s) specified
%   by DIM.  For example, SIZE(X,1) returns the size of the first
%   dimension of the tensor.
%
%   See also TENMAT, TENMAT/SIZE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isempty(a.data)
    sz = [];
elseif exist('idx', 'var')
    sz = a.tsize(idx);
else
    sz = a.tsize;
end

return;
