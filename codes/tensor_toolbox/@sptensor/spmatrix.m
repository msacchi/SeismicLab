function s = spmatrix(a)
%SPMATRIX Converts a two-way sparse tensor to sparse matrix.
%
%   SPMATRIX(X) converts a sparse tensor to a sparse matrix. The sparse
%   tensor must be two-dimensional.
%
%   See also SPTENSOR, SPTENSOR/RESHAPE, SPTENMAT
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if ndims(a) ~= 2
    error('Sparse tensor must be two dimensional.');
end


if isempty(a.subs)
    s = sparse(a.size(1), a.size(2));
else
    s = sparse(a.subs(:,1), a.subs(:,2), a.vals, a.size(1), a.size(2));
end
