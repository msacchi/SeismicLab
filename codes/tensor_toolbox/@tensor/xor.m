function Z = xor(X,Y)
%XOR Logical EXCLUSIVE OR for tensors.
%
%   See also TENSOR.
% 
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



Z = tenfun(@xor,X,Y);
