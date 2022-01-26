function Z = ldivide(X,Y)
%LDIVIDE Left array divide for tensor.
%
%   LDIVIDE(A,B) is called for the syntax 'A .\ B' when A or B is a tensor.
%   A and B must have the same size, unless one is a scalar.  
%
%   Examples
%   X = tenrand([4 3 2],5);
%   X .\ 3
%   X .\ X
%
%   See also TENSOR, TENSOR/RDIVIDE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



Z = tenfun(@ldivide,X,Y);
