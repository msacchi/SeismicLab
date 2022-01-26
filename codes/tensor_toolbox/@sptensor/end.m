function e = end(X,k,n)
%END Last index of indexing expression for sparse tensor.
%
%   The expression X(end,:,:) will call END(X,1,3) to determine
%   the value of the first index.
%
%   See also SPTENSOR, SPTENSOR/SUBSREF, END.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


switch n
    case 1 %linear indexing
        e = prod(X.size);
    case ndims(X) %subscript indexing
        e = X.size(k);
    otherwise
        error('Invalid subscripting');
end
    
