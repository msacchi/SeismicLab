function c = ttv(a,v,dims)
%TTV Tensor times vector for ttensor.
%
%   Y = TTV(X,A,N) computes the product of ttensor X with a
%   (column) vector A.  The integer N specifies the dimension in X
%   along which A is multiplied.  If size(A) = [I,1], then X must have
%   size(X,N) = I.  Note that ndims(Y) = ndims(X) - 1 because the N-th
%   dimension is removed.
%
%   Y = TTV(X,{A1,A2,...}) computes the product of ttensor X with a
%   sequence of vectors in the cell array.  The products are computed
%   sequentially along all dimensions (or modes) of X. The cell array
%   contains ndims(X) vectors.
%
%   Y = TTV(X,{A1,A2,...},DIMS) computes the sequence tensor-vector
%   products along the dimensions specified by a vector DIMS.
%
%   Examples
%   X = ttensor(tensor(rand(2,3,2)), rand(3,2), rand(2,3), rand(2,2));
%   ttv(X, [1:3]', 1)
%   ttv(X, {[1:3]', [1:2]', [1:2]'})
%   ttv(X, {[1:3]', [1:2]'}, [1 3])
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','multiply_doc.html')))">Documentation page for multiplying tensors</a>
%
%   See also TENSOR/TTV, TTENSOR, TTENSOR/TTM.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

%%%%%%%%%%%%%%%%%%%%%%
%%% ERROR CHECKING %%%
%%%%%%%%%%%%%%%%%%%%%%

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

% Check that each multiplicand is the right size.
for i = 1:numel(dims)
    if ~isequal(size(v{vidx(i)}),[size(a,dims(i)) 1])
        error('Multiplicand is wrong size');
    end
end

% Figure out which dimensions will be left when we're done
remdims = setdiff(1:ndims(a),dims);

% Create w to be multiplied with a.core
w = cell(ndims(a),1);
for i = 1:numel(dims)
    w{dims(i)} = a.u{dims(i)}' * v{vidx(i)};
end

% Create new core
newcore = ttv(a.core,w,dims);

% Create final result
if isempty(remdims)
    c = newcore;
else
    c = ttensor(newcore,a.u{remdims});
end
