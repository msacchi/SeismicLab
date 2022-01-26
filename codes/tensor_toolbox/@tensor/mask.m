function vals = mask(X,W)
%MASK Extract values as specified by a mask tensor.
%
%   V = MASK(X,W) extracts the values in X that correspond to nonzero
%   values in the mask tensor W.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

% Error check
if any(size(W) > size(X))
    error('Mask cannot be bigger than the data tensor')
end

% Extract locations of nonzeros in W
wsubs = find(W);

% Extract values from X
idx = tt_sub2ind(X.size,wsubs);
vals = X.data(idx);
