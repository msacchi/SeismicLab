function Z = ldivide(X,Y)
%LDIVIDE Left array divide (.\) for symmetric tensors.
%
%   LDIVIDE(A,B) is called for the syntax 'A .\ B' when A or B is a symmetric 
%   tensor. A and B must have the same size, unless one is a scalar.  
%
%   Examples
%   X = symtenrand([4 4 4]);
%   X .\ 3
%   X .\ X
%
%   See also SYMTENSOR, SYMTENSOR/RDIVIDE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



Z = tenfun(@ldivide,X,Y);
