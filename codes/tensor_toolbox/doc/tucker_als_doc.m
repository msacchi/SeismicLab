%% Alternating least squares for Tucker model 
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="tucker.html">Tucker Decompositions</a> 
% &#62;&#62; <a href="tucker_als_doc.html">Tucker-ALS</a>
% </p>
% </html>
%
% The function |tucker_als| computes the best rank(R1,R2,..,Rn)
% approximation of tensor X, according to the specified dimensions in
% vector R.  The input X can be a tensor, sptensor, ktensor, or
% ttensor.  The result returned in T is a ttensor.
%
% The method is originally from Tucker (1966) and later revisited in 
% De Lathauwer et al. (2000).
%
% * L. R. Tucker, Some mathematical notes on three-mode factor analysis,
%   Psychometrika, 31:279-311, 1966, http://dx.doi.org/10.1007/BF02289464
% * L. De Lathauwer, B. De Moor, J. Vandewalle,
%   On the best rank-1 and rank-(R_1, R_2, R_N) approximation of
%   higher-order tensors,
%   SIAM J. Matrix Analysis and Applications, 21:1324-1342, 2000, 
%   http://doi.org/10.1137/S0895479898346995  
%
% Note: Oftentimes it's better to use |hosvd| instead.

%% Create a data tensor of size [5 4 3]
rng('default'); rng(0); %<-- Set seed for reproducibility
X = sptenrand([5 4 3], 10)
%% Create a [2 2 2] approximation
T = tucker_als(X,2)        %<-- best rank(2,2,2) approximation 
%% Create a [2 2 1] approximation
T = tucker_als(X,[2 2 1])  %<-- best rank(2,2,1) approximation 
%% Use a different ordering of the dimensions
T = tucker_als(X,2,struct('dimorder',[3 2 1]))
%% Use the n-vecs initialization method
% This initialization is more expensive but generally works very well.
T = tucker_als(X,2,struct('dimorder',[3 2 1],'init','eigs'))
%% Specify the initial guess manually
U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of T
T = tucker_als(X,2,struct('dimorder',[3 2 1],'init',{U0}))
