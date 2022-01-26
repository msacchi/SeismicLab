function a = double(s)
%DOUBLE Converts a sparse tensor to a dense multidimensional array.
%
%   See also SPTENSOR, SPTENSOR/FULL.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



a = zeros([size(s) 1 1]);
if nnz(s) > 0
    a(tt_sub2ind(size(s),s.subs)) = s.vals;
end
