function new_X = extract(X,idx)
%EXTRACT Creates a new ktensor with only the specified components.
%
%   Y = EXTRACT(X,S) selected the subset of components in X as defined by
%   S. It should be the case that S is a subset of [1,...,NCOMPONENTS(X)].
%
%   See also KTENSOR, NCOMPONENTS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%% Set-up
N = ndims(X);
%% Extract
new_lambda = X.lambda(idx);
new_U = cell(N,1);
for i = 1 : N
    new_U{i} = X.u{i}(:,idx);
end
new_X = ktensor(new_lambda, new_U);

