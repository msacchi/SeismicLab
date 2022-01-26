function A = double(T)
%DOUBLE Convert a sptenmat to a sparse matrix.
%
%   A = double(T) converts T stored as a SPTENMAT to a sparse matrix.
%
%   See also SPTENMAT.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



m = prod(T.tsize(T.rdims));
n = prod(T.tsize(T.cdims));
if isempty(T.subs)
    A = sparse(m,n);
else
    A = sparse(T.subs(:,1), T.subs(:,2), T.vals, m, n);
end
