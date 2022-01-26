function c = ttv(a,v,dims)
%TTV Tensor times vector for sumtensor.
%
%   Y = TTV(X,A,N) computes the product of sumtensor X with a
%   (column) vector A.  The integer N specifies the dimension in X
%   along which A is multiplied.  If size(A) = [I,1], then X must have
%   size(X,N) = I.  Note that ndims(Y) = ndims(X) - 1 because the N-th
%   dimension is removed.
%
%   Y = TTV(X,{A1,A2,...}) computes the product of sumtensor X with a
%   sequence of vectors in the cell array.  The products are computed
%   sequentially along all dimensions (or modes) of X. The cell array
%   contains ndims(X) vectors.
%
%   Y = TTV(X,{A1,A2,...},DIMS) computes the sequence tensor-vector
%   products along the dimensions specified by DIMS.
%
%   Examples
%   T1 = tensor(rand(3,3,3));
%   T2 = sptensor([1 1 1; 3 1 2; 1 1 3], 1, [3,3,3]);
%   T = sumtensor(T1, T2); %<--Declaring a sumtensor
%
%   ttv(T, [1 1 1]', 3) %<--Multiply ones vector along mode 3
%   ttv(T, {[1 1 1]', [1 1 1]', [1 1 1]'}) %<--ones vector along all modes
%   ttv(T, {[1 1 1]', [1 1 1]'}, [1 3]) %<--ones vector along modes 1 and 3
%
%   See also TENSOR/TTV, SUMTENSOR, SUMTENSOR/TTM.
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

tmp = cell(length(a.part),1);
for i = 1:length(a.part)
    tmp{i} = ttv(a.part{i},v,dims);
end

if isscalar(tmp{1})
    c = sum(cell2mat(tmp));
else
    c = sumtensor(tmp{:});
end


