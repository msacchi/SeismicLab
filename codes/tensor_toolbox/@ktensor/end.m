function e = end(X,k,n)
%END Last index of indexing expression for ktensor.
%
%   The expression X(end,:,:) will call END(X,1,3) to determine
%   the value of the first index.
%
%   See also KTENSOR, KTENSOR/SUBSREF, END.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%TODO (after 2.0 release): Resolve ambiguity w.r.t X{end}and X(end,1,1)
%for 1st-order tensors.

if n > ndims(X)
  error('Subscript out of range.');
end

if (n ~= 1)
  e = size(X,k);
else
  e = ndims(X);
end
