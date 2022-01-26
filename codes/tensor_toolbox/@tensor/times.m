function Z = times(X,Y)
%TIMES Array multiplication for tensors.
%
%   TIMES(A,B) is called for the syntax 'A .* B' when A or B is a 
%   tensor. A and B must have the same size, unless one is a scalar. A
%   scalar can be multiplied by a tensor of any size.
%
%   See also TENSOR.
% 
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isa(Y,'ktensor') || isa(Y,'ttensor') || isa(Y,'sptensor')
    Y = full(Y);
end

Z = tenfun(@times,X,Y);
