function Z = and(X,Y)
%AND Logical AND (&) for symmetric tensors.
%
%   See also SYMTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


  
Z = tenfun(@and,X,Y);
