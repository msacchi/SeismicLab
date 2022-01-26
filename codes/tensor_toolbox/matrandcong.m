function X = matrandcong(m,n,gamma)
%MATRANDCONG Create a random matrix with a fixed congruence.
%
%   X = MATRANDCONG(M,N,GAMMA) creates a matrix X of size M x N such
%   that each column of X has norm 1 and any two columns of X have an inner
%   product equal to GAMMA.
%
%   Based on code from Evrim Acar and the paper G. Tomasi and R. Bro, A
%   comparison of algorithms for fitting the PARAFAC model, Computational
%   Statistics & Data Analysis, 50: 1700-1734, 2006.
%
%   See also MATRANDORTH, MATRANDNORM, CREATE_PROBLEM, CREATE_GUESS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

CG = gamma * ones(n,n) + (1-gamma) * eye(n);
CGR = chol(CG);
X = randn(m,n);
[Q,~] = qr(X,0);
X = Q * CGR;
