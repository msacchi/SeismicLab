function [tf,diffs] = issymmetric(X)
%ISSYMMETRIC Verify that a ktensor X is symmetric in all modes.
%
%   TF = ISSYMMETRIC(X) returns true if X is exactly symmetric for every
%   permutation.
%
%   [TF,DIFFS] = ISSYMMETRIC(X) also returns the matrix of the norm of the
%   differences between the normalized factor matrices.
%
%   See also SYMMETRIZE.
%  
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%T. Kolda, June 2014.

n = ndims(X);
sz = size(X);
diffs = zeros(n,n);

for i = 1:n
    for j = i+1:n
        if ~isequal(size(X.u{i}), size(X.u{j}))
            diffs(i,j) = Inf;            
        elseif isequal(X.u{i},X.u{j})
            diffs(i,j) = 0;
        else
            diffs(i,j) = norm(X.u{i} - X.u{j});
        end
    end
end

tf = all(diffs(:) == 0);
