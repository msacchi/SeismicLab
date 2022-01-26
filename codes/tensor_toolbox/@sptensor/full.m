function B = full(A)
%FULL Convert a sparse tensor to a (dense) tensor.
%
%   B = FULL(A) converts a sptensor A to a (dense) tensor B.
%
%   See also SPTENSOR, TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



% Extract the order and size of A
siz = size(A);

% Handle the completely empty (no size) case
if isempty(siz)
    B = tensor;
    return;
end

% Create a dense zero tensor B that is the same size as A
B = tensor(zeros([siz,1,1]),siz);

if isempty(A.subs)
    return;
end

% Extract the linear indices of entries in A
idx = tt_sub2ind(siz,A.subs);

% Copy the values of A into B using linear indices
B(idx) = A.vals;

