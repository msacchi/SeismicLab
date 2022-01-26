function Z = plus( X, Y )
%PLUS Binary addition (+) for symtensors.
%
%   PLUS(A,B) is called for the syntax 'X + Y' when X and Y are symtensors.
%   A and B must be the same size, unless one is a scalar. A scalar can be
%   added to a symtensor of any size.
%
%   See also SYMTENSOR, TENSOR/PLUS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


Z = tenfun(@plus,X,Y);


    


