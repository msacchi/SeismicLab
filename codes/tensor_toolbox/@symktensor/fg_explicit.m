function [f,g] = fg_explicit(M,X,normXsqr)
%FG_EXPLICIT Function and gradient of F(M)=||X-M||^2 for symktensor.
%
%   [F,G] = FG_EXPLICIT(M,X,NORMXSQR) requires that M be a symktensor and X
%   be a (dense) tensor that is symmetric. It computes the function (F) and
%   corresponding gradient (G) of 
%
%      F(M) = ||X-M||^2 = ||X||^2 - 2 <M,X> + ||M||^2
%
%   where the value of ||X||^2 is provided as NORMXSQR. Because this term
%   does not depend on M, it can be ignored (replaced by zero) in the
%   context of optimization. The gradient G is returned in vectorized form.
%   If M is a d-way, n-dimensional, rank-r symktensor, then G is a vector
%   of length (n+1)*r.    
%
%   See also symktensor/fg, symktensor/fg_implicit.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

%Code by Samantha Sherman, Summer 2018.

lambda = M.lambda;
A = M.u;

d = ndims(X);
r = size(A,2);

Y = zeros(size(A));
for j = 1:r
    Y(:,j) = ttsv(X, A(:,j), -1);    % X * a^(d-1) 
end

w = bsxfun(@dot,A,Y)';

B = A'*A;
C = B.^(d-1); %(B^TB)^(d-1)
u = (C.*B)*lambda;

f = normXsqr + dot(lambda,u) - 2 * dot(lambda,w);

G.lambda = -2 * (w - u);

G.A = -2 * d * (Y - A * diag(lambda) * C) * diag(lambda);

g = [G.lambda; reshape(G.A,[],1)];


