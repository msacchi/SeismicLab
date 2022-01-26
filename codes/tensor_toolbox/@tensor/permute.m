function T = permute(T,order)
%PERMUTE Permute tensor dimensions.
%
%   B = PERMUTE(A,ORDER) rearranges the dimensions of A so that they
%   are in the order specified by the vector ORDER. The result has the
%   same values of A, but the order of the subscripts needed to access
%   any particular element are rearranged as specified by ORDER.
%
%   Examples
%   T = tensor(rand(3,2,4));
%   permute(T,[1 3 2])
%
%   See also TENSOR, TENSOR/SIZE, TENSOR/NDIMS, PERMUTE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if ndims(T) ~= numel(order)
  error('Invalid permutation order');
end

% Check for special case of permuting an order-1 object (which has
% no effect but confuses MATLAB's permute command which doesn't
% think that there is such a thing as a 1D-array).
if isequal(order,1)
    return;
end

% Check for special case of empty object (which has
% no effect but confuses MATLAB's permute command which doesn't
% think that there is such a thing as an empty array).
if isempty(order)
    return;
end

% Note that permute does error checking on order, so we don't worry
% about it. 
T.data = permute(T.data,order);
T.size = T.size(order);

return;
