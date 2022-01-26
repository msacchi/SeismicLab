function f = f_implicit(M,X,normXsqr)
%F_IMPLICIT Function F(M)=||X-M||^2 for two symktensors.
%
%   F = F_IMPLICIT(M,X,NORMXSQR) requires *both* M and X be symktensor
%   objects.  It computes the function value
%
%      F(M) = ||X-M||^2 = ||X||^2 - 2 <M,X> + ||M||^2
%
%   where the value of ||X||^2 is provided as NORMXSQR. Because this term
%   does not depend on M, it can be ignored (replaced by zero) in the
%   context of optimization. 
%
%   See also CP_ISYM, SYMKTENSOR/FG_IMPLICIT, SYMKTENSOR/G_IMPLICIT.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

%% extract data

d = ndims(M);
V = X.u;
eta = X.lambda;

A = M.u;
lambda = M.lambda;

%% do computations
M1 = (V'*A).^(d-1);
M2 = bsxfun(@times,V',eta);
Y = M2' * M1;
w = bsxfun(@dot,A,Y)';

B = A'*A;
u = (B.^d)*lambda;

f = normXsqr + dot(lambda,u) - 2 * dot(lambda,w);

