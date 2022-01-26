function c = ttv(a,v,dims)
%TTV Sparse tensor times vector.
%
%   Y = TTV(X,V,N) computes the product of a sparse tensor X with a
%   (column) vector V.  The integer N specifies the dimension in X
%   along which V is multiplied.  If size(V) = [I,1], then X must have
%   size(X,N) = I.  Note that ndims(Y) = ndims(X) - 1 because the N-th
%   dimension is removed. 
%
%   Y = TTV(X,U) computes the product of a sparse tensor X with a
%   sequence of vectors in the cell array U.  The products are
%   computed sequentially along all dimensions (or modes) of X. The
%   cell array U contains ndims(X) vectors.
%
%   Y = TTV(X,U,DIMS) computes the sequence tensor-vector products
%   along the dimensions specified by DIMS.
%
%   In all cases, the result Y is a sparse tensor if it has 50% or
%   fewer nonzeros; otherwise the result is returned as a dense
%   tensor.
%
%   See also SPTENSOR, SPTENSOR/TTM, TENSOR, TENSOR/TTV.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



% Check the number of arguments
if (nargin < 2)
    error('TTV requires at least two arguments.');
end

% Check for 3rd argument
if ~exist('dims','var')
    dims = [];
end

% Check that 2nd argument is cell array. If not, recall with v as a
% cell array with one element.
if ~iscell(v)
    c = ttv(a,{v},dims);
    return;
end

% Get sorted dims and index for multiplicands
[dims,vidx] = tt_dimscheck(dims,ndims(a),numel(v));       
remdims = setdiff(1:ndims(a),dims);

% Check that each multiplicand is the right size.
for i = 1:numel(dims)
    if ~isequal(size(v{vidx(i)}),[size(a,dims(i)) 1])
        error('Multiplicand is wrong size');
    end
end

% Multiply each value by the appropriate elements of the
% appropriate vector
newvals = a.vals;
subs = a.subs;
if isempty(subs) %There are no nonzero terms in a
    newsubs = [];
else
    for n = 1:length(dims)
        idx = subs(:,dims(n)); % extract indices for dimension n
        w = v{vidx(n)};        % extract nth vector
        bigw = w(idx);         % stretch out the vector
        newvals = newvals .* bigw;
    end
    newsubs = subs(:,remdims);
end
% Case 0: If all dimensions were used, then just return the sum
if isempty(remdims)
    c = sum(newvals);
    return;
end

% Otherwise, figure new subscripts and accumulate the results.
newsiz = a.size(remdims);

% Case I: Result is a vector
if numel(remdims) == 1
    c = accumarray(newsubs,newvals,[newsiz 1]);
    if nnz(c) <= 0.5 * newsiz
        c = sptensor((1:newsiz)',c,newsiz);
    else
        c = tensor(c,newsiz);
    end
    return;
end

% Case II: Result is a multiway array
c = sptensor(newsubs, newvals, newsiz);

% Convert to a dense tensor if more than 50% of the result is nonzero.
if nnz(c) > 0.5 * prod(c.size)
    c = tensor(c);
end

return;
