function Z = rdivide(X,Y)
%RDIVIDE Right array divide (./) for tensors.
%
%   RDIVIDE(A,B) is called for the syntax 'A ./ B' when A or B is a tensor.
%   A and B must have the same size, unless one is a scalar.  
%
%   Examples
%   X = tenrand([4 4 4]);
%   X ./ 3
%   X ./ X
%
%   See also SYMTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


Z = tenfun(@rdivide,X,Y);
