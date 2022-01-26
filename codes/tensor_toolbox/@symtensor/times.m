function Z = times(X,Y)
%TIMES Array multiplication (.*) for symmetric tensors.
%
%   TIMES(A,B) is called for the syntax 'A .* B' when A or B is a 
%   symmetric tensor. A and B must have the same size, unless one is a 
%   scalar. A scalar can be multiplied by a symmetric tensor of any size.
%
%   See also SYMTENSOR.
% 
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


Z = tenfun(@times,X,Y);
