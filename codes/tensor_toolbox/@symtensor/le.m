function Z = le(X,Y)
%LE Less than or equal (<=) for symmetric tensors.
%
%   See also SYMTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



Z = tenfun(@le,X,Y);
