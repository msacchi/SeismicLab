function sz = size(a,idx)
%SIZE Size of tenmat.
%
%   D = SIZE(X) returns the two-element row vector D = [M N]
%   containing the number of rows and columns in the matrix.
% 
%   M = SIZE(X,DIM) returns the length of the dimension specified by
%   the scalar DIM.  For example, SIZE(X,1) returns the number of
%   rows.
%
%   See also TENMAT, TENMAT/TSIZE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isempty(a.data)
    sz = [];
elseif exist('idx', 'var')
    sz = size(a.data, idx);
else
    sz = size(a.data);
end
