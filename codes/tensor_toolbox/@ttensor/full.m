function X = full(T)
%FULL Convert a ttensor to a (dense) tensor.
%
%   X = FULL(T) converts ttensor T to (dense) tensor X.
%
%   See also TTENSOR, TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


% Preallocate to ensure there is enough space
X = tenzeros(size(T));

% Now do the calculation 
X = ttm(T.core,T.u);

% Make sure that X is a dense tensor (small chance it could be a sparse
% tensor).
X = tensor(X);

return;
