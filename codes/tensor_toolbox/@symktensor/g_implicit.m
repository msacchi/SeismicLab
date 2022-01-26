function g = g_implicit(M,X)
%G_IMPLICIT Gradient of F(M)=||X-M||^2 for two symktensors.
%
%   G = G_IMPLICIT(M,X,NORMXSQR) requires *both* M and X be symktensor
%   objects and computes the gradient (G) of the function
%
%      F(M) = ||X-M||^2 = ||X||^2 - 2 <M,X> + ||M||^2.
%
%   The gradient G is returned in vectorized form. If M is a d-way,
%   n-dimensional, rank-r symktensor, then G is a vector of length (n+1)*r.    
%
%   See also cp_isym, symktensor/fg_implicit, symktensor/f_implicit.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

%Code by Samantha Sherman, Summer 2018.

% extract data

d = ndims(X);
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
C = B.^(d-1); %(B^TB)^(d-1)
u = (C.*B)*lambda;

G.lambda = -2 * (w - u);

G.A = -2 * d * (Y - A * diag(lambda) * C) * diag(lambda);

g = [G.lambda; reshape(G.A,[],1)];


