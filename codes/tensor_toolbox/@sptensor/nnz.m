function a = nnz(t)
%NNZ Number of nonzeros in sparse tensor.
%
%   NNZ(T) is the number of nonzero elements in T.
%
%   See also SPTENSOR, SPTENSOR/FIND.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isempty(t.subs)
    a = 0;
else
    a = size(t.subs,1);
end
