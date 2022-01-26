function y = contract(x,i,j)
%CONTRACT Contract sparse tensor along two dimensions (array trace).
%
%   Y = CONTRACT(X,I,J) contracts the entries of X along dimensions I
%   and J. Contraction is a generalization of matrix trace. In other
%   words, the trace is performed along the two-dimensional slices
%   defined by dimensions I and J. It is possible to implement tensor
%   multiplication as an outer product followed by a contraction.
%
%   Examples
%   X = sptenrand([4 3 2],10); Y = sptenrand([3 2 4],10);
%   Z1 = ttt(X,Y,1,3); %<-- Normal tensor multiplication
%   Z2 = contract(ttt(X,Y),1,6); %<-- Outer product + contract
%   norm(Z1-Z2) %<-- Should be zero
%
%   See also SPTENSOR, SPTENSOR/TTT.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



% Error checking
if x.size(i) ~= x.size(j)
    error('Must contract along equally sized dimensions');
end

% Error checking
if i == j
    error('Must contract along two different dimensions');
end

% Easy case - returns a scalar
if ndims(x) == 2
    tfidx = (x.subs(:,1) == x.subs(:,2)); % find diagonal entries
    y = sum(x.vals(tfidx));
    return;
end

% Remaining dimensions after contract
remdims = setdiff(1:ndims(x),[i j]);

% Find index of values on diagonal
indx = find(x.subs(:,i) == x.subs(:,j));

% Let the constructor sum up the entries
y = sptensor(x.subs(indx,remdims),x.vals(indx),x.size(remdims));

% Check if result should be dense
if nnz(y) > 0.5 * prod(y.size)
    % Final result is a *dense* tensor
    y = tensor(y);
end
