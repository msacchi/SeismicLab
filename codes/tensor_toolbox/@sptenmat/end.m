function e = end(X,k,n)
%END Last index of indexing expression for sptenmat.
%
%   The expression X(end,:) will call END(X,1,2) to determine
%   the value of the first index.
%
%   See also SPTENMAT, SPTENMAT/SUBSREF, END.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if n > 2
  error('Subscript out of range.');
end
e = size(X,k);
