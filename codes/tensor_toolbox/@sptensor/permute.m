function t = permute(t,order)
%PERMUTE Rearrange the dimensions of a sparse tensor.
%
%   B = PERMUTE(A,ORDER) rearranges the dimensions of A so that they
%   are in the order specified by the vector ORDER. The result has the
%   same values of A, but the order of the subscripts needed to access
%   any particular element are rearranged as specified by ORDER.
%
%   See also SPTENSOR, PERMUTE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



% Error checking
if (ndims(order) ~= 2) || (size(order,1) ~= 1) 
    error('ORDER must be a row vector');
end
   
% Check that the permuation is valid
if ~isequal(sort(order),1:ndims(t))
    error('Invalid permutation.');
end

% Do the permutation
if ~isempty(t.subs)
    t.subs = t.subs(:,order);
end
t.size = t.size(order);
