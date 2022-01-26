function X = matrandnorm(varargin)
%MATRANDNORM Normalizes columns of X so that each is unit 2-norm.
%
%   X = MATRANDNORM(M,N) creates a random M x N matrix with randomly using
%   normally distributed enries and then rescales the columsn so that each
%   has a unit 2-norm.
%
%   X = MATRANDNORM(X) rescales the columns of X so that each
%   column has a unit 2-norm. 
%
%   Examples
%      X = MATRANDNORM(rand(5,5));
%      X = MATRANDNORM(3,2);
%      X = MATRANDNORM(ones(4));
% 
%   See also MATRANDORTH, MATRANDNORM, CREATE_PROBLEM, CREATE_GUESS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

if nargin == 2
    X = randn(varargin{1}, varargin{2});
else
    X = varargin{1};
end

norms = sqrt(sum(X.^2,1));
X = bsxfun(@rdivide,X,norms);

