function Z = or(X,Y)
%OR Logical OR (|) for symmetric tensors.
%
%   See also SYMTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



Z = tenfun(@or,X,Y);
