function B = full(A)
%FULL Convert a sptenmat to a (dense) tenmat.
%
%   B = FULL(A) converts a sptenmat A to a (dense) tenmat B.
%
%   See also SPTENMAT, TENMAT.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



% Extract the order and size of A
siz = size(A);

% Create a dense zero tensor B that is the same size as A
B = tenmat(zeros([siz,1,1]), A.rdims, A.cdims, A.tsize);

% Extract the linear indices of entries in A
idx = tt_sub2ind(siz,A.subs);

% Copy the values of A into B using linear indices
B(idx) = A.vals;
