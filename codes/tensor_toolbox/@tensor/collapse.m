function Y = collapse(X,dims,fun)
%COLLAPSE Collapse tensor along specified dimensions.
%
%   Y = COLLAPSE(X,DIMS) sums the entries of X along all dimensions
%   specified in DIMS. If DIMS is negative, then X is summed across
%   all dimensions *not* specified by -DIMS.
%
%   Y = COLLAPSE(X) is shorthand for S = COLLAPSE(X,1:ndims(X)).
%
%   Y = COLLAPSE(X,DIMS,FUN) accumulates the entries of T using the
%   accumulation function @FUN.
%
%   Examples
%   X = tenrand([4 4 4]);
%   Y = collapse(X,[2 3]) %<-- sum of entries in each mode-1 slice
%   Y = collapse(X,[1 2],@max) %<-- max entry in each mode-3 slice
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','collapse_scale_doc.html')))">Documentation page for collapsing and scaling tensors</a>
%
%   See also TENSOR, TENSOR/SCALE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isempty(X.data)
    Y = [];
    return;
end

if ~exist('dims', 'var')
    dims = 1:ndims(X);
end

if isempty(dims)
    Y = X;
    return;
end

if ~exist('fun', 'var')
    fun = @sum;
end

dims = tt_dimscheck(dims,ndims(X));
remdims = setdiff(1:ndims(X),dims);

% Check for the case where we accumulate over *all* dimensions
if isempty(remdims)
    Y = fun(X.data(:));
    return;
end

% Calculate the size of the result
newsiz = size(X,remdims);

% Convert to a matrix where each row is going to be collapsed
A = double(tenmat(X,remdims,dims));

% Apply the collapse function
B = zeros(size(A,1),1);
for i = 1:size(A,1)
    B(i) = fun(A(i,:));
end

% Form and return the final result
Y = tensor(tenmat(B,1:numel(remdims),[],newsiz));




