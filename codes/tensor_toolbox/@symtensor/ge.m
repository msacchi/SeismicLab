function Z = ge(X,Y)
%GE Greater than or equal (>=) for tensors.
%
%   See also SYMTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



Z = tenfun(@ge,X,Y);
