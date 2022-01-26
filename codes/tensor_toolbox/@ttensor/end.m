function e = end(X,k,n)
%END Last index of indexing expression for ttensor.
%
%   The expression X(end,:,:) will call END(X,1,3) to determine
%   the value of the first index.
%
%   See also TTENSOR, TTENSOR/SUBSREF, END.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


% Note that this only works with {} because () is not supported by
% subsref.
e = ndims(X);
