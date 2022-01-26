function e = end(X,k,n)
%END Last index of indexing expression for tensor.
%
%   The expression X(end,:,:) will call END(X,1,3) to determine
%   the value of the first index.
%
%   See also TENSOR, TENSOR/SUBSREF, END.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if n > ndims(X)
  error('Subscript out of range.');
end
if n > 1 %For subscripted indexing
    e = X.size(k); %For subscripted indexing
else %Linear indexing, or X is a vector
    e = prod(size(X)); %if X is a vector, this equals X.size(1) so works
end
    
