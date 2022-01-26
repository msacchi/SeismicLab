function n = ncomponents(t)
%NCOMPONENTS Number of components for a symktensor.
%
%   NCOMPONENTS(T) returns the number of components in the symktensor T.
%   This is the size of the lambda vector, equivalently the number of
%   columns in the factor matrix.
%
%   S = symktensor(3, symtensor(@rand,4,3)); %<--Random symktensor
%   ncomponents(S) %<--Returns 3
%
%   See also SYMKTENSOR
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



n = length(t.lambda);
