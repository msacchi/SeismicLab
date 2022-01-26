function Z = minus(X,Y)
%MINUS Binary subtraction (-) for tensors.
%
%   MINUS(A,B) is called for the syntax 'A - B' when A or B is a tensor. A
%   and B must have the same size, unless one is a scalar. A scalar can be
%   subtracted from a tensor of any size. 
%
%   See also TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



Z = tenfun(@minus,X,Y);
