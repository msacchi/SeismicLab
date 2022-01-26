function X = arrange(X,perm)
%ARRANGE Arranges the rank-1 components of a symktensor.
%
%   ARRANGE(X) normalizes the columns of the factor matrices and then sorts
%   the components by magnitude, greatest to least.
%
%   ARRANGE(X,P) rearranges the components of X according to the
%   permutation P. P should be a permutation of 1 to NCOMPONENTS(X). 
%
%   See also SYMKTENSOR, NCOMPONENTS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%% Just rearrange and return if second argument is a permutation
if exist('perm','var') && isvector(perm)
    X.lambda = X.lambda(perm);
    X.u = X.u(:,perm);
    return;
end

%% Ensure that matrices are normalized
X = normalize(X);

%% Sort
[~, idx] = sort(abs(X.lambda), 1, 'descend');
X.lambda = X.lambda(idx);
X.u = X.u(:,idx);


