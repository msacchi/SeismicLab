function Z = xor(X,Y)
%XOR Logical EXCLUSIVE OR for symmetric tensors.
%
%   See also SYMTENSOR.
% 
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



Z = tenfun(@xor,X,Y);
