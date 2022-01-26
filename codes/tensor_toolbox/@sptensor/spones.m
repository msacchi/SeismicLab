function X = spones(X)
%SPONES Replace nonzero sparse tensor elements with ones.
%
%   Y = SPONES(X) generates a tensor with the same sparsity structure as X,
%   but with ones in the nonzero positions.
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


X.vals = ones(size(X.vals));
