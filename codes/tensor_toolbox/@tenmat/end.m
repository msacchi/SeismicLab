function e = end(X,k,n)
%END Last index of indexing expression for tenmat.
%
%   The expression X(end,:) will call END(X,1,2) to determine
%   the value of the first index.
%
%   See also TENMAT, TENMAT/SUBSREF, END.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if n > ndims(X)
  error('Subscript out of range.');
end
e = size(X.data,k);
