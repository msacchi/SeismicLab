function Z = plus(X,Y)
%PLUS Binary addition (+) for tensors. 
%
%   PLUS(A,B) is called for the syntax 'A + B' when A or B is a tensor. A
%   and B must have the same size, unless one is a scalar. A scalar can be
%   added to a tensor of any size. 
%
%   See also TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if isa(Y,'sumtensor') %If the 2nd component is a sumtensor, treat as such
    Z=plus(Y,X)
else
    Z = tenfun(@plus,X,Y);
end
