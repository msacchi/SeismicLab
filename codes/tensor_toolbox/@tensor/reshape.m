function t = reshape(t,siz)
%RESHAPE Change tensor size.
%
%   RESHAPE(X,SIZ) returns the tensor whose elements
%   have been reshaped to the appropriate size.
%
%   Examples
%   X = tensor(rand(2,3,4));
%   reshape(X,[4,3,2])
%
%   See also TENSOR, TENSOR/SQUEEZE, TENSOR/PERMUTE.
%   
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if prod(t.size) ~= prod(siz)
    error('Number of elements cannot change');
end

t.data = reshape(t.data,siz);
t.size = siz;
