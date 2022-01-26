%% Shifted Power Method for Generalized Tensor Eigenproblem
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="eigen.html">Tensor Eigenproblem</a> 
% &#62;&#62; <a href="sshopm_doc.html">GEAP</a>
% </p>
% </html>
%
% The method is described in the following paper:
%
%   * T. G. Kolda, J. R. Mayo, An Adaptive Shifted Power Method for
%     Computing Generalized Tensor Eigenpairs, SIAM J. Matrix Analysis and
%     Applications, 35:1563-1582, 2014, http://dx.doi.org/0.1137/140951758    
%

%%
rng('default'); %<- Setup for reproducibility

%% Data tensor for Example 5.1 in Kolda and Mayo (2014)
% Tensor taken from Example 1 in E. Kofidis and P. A. Regalia, On the best
% rank-1 approximation of higher-order supersymmetric tensors, SIAM J.
% Matrix Anal. Appl., 23:863-884, 2002,
% http://dx.doi.org/10.1137/S0895479801387413. 

% Unique values in lexiographical order
Avals = [0.2883 -0.0031 0.1973 -0.2485 -0.2939 0.3847 0.2972 0.1862 ...
    0.0919 -0.3619 0.1241 -0.3420 0.2127 0.2727 -0.3054 ];

% Create standard tensor object directly
% because that's what EIG_GEAP expects.
A = full(symtensor(Avals, 4,3)); 

% Display the tensor as a symmetric tensor
disp(symtensor(A),'A')

% Save order and dimension
m = ndims(Asym); 
n = size(Asym,1); 

%% Create corresponding "identity" tensor 
% Create an identity tensor of order m and dimension n.
% Here we use a somewhat convoluted construction method, but it aligns with
% the mathematical derivation.
if mod(m,2) ~= 0
    error('Identity tensor on exists for even order');
end
Bsym = symtensor(@zeros,m,n);
uniqidx = indices(Bsym);
Bvals = zeros(size(uniqidx,1),1);
for i = 1:size(uniqidx,1)
    pidx = perms(uniqidx(i,:));
    pidxodd = pidx(:,1:2:m-1);
    pidxeven = pidx(:,2:2:m);
    pidxresult = pidxodd == pidxeven;
    Bvals(i) = sum(all(pidxresult,2))/factorial(m);
end
Bsym = symtensor(Bvals, m,n);
B = full(Bsym);
disp(Bsym,'B')

%% Demonstrate identity property
% A _tensor identity_, E,  satisfies the following mathematical property for any
% n-dimensional vector x. 
% 
% $$\|x\|^2 x = E x^{\otimes 2}$$
%
% Note that it is not scale invariant. When we test this property, it
% should be close to machine precision.
% 

for i = 1:10
    x = rand(n,1);
    lhs = norm(x)^(m-2) * x;
    rhs = ttsv(B,x,-1);
    fprintf('Identity property error: %g\n', norm(lhs-rhs));
end

%% Call |eig_geap| to find eigenpair
% Default is to find a maxima, i.e., a convex solution.
info = eig_geap(A, B, 'MaxIts', 100, 'Display',1);

%% Demonstrate generalized eigenpair property
% 
% $$A x^{\otimes 2} = \lambda B x^{\otimes 2}$$
% 

x = info.x;
lambda = info.lambda;
lhs = ttsv(A,x,-1);
rhs = lambda * ttsv(B,x,-1);
fprintf('Generalized eigenpair identity norm: %g\n', norm(lhs-rhs));

%% Call |eig_geap| to find another eigenpair
% Find a minima (by specifying 'Concave' = true)
info = eig_geap(A, B, 'MaxIts', 100, 'Concave', true, 'Display',1);

%% 
% Test generalized eigenvector properties
x = info.x;
lambda = info.lambda;
lhs = ttsv(A,x,-1);
rhs = lambda * ttsv(B,x,-1);
fprintf('Generalized eigenpair identity norm: %g\n', norm(lhs-rhs));

%% Reproduce tensors from Example 5.3 from Kolda and Mayo (2014)

Avals = [0.4982 -0.0582 -1.1719 0.2236 ...
    -0.0171 0.4597 0.4880 0.1852 ...
    -0.4087 0.7639 0.0000 -0.6162 ...
    0.1519 0.7631 2.6311];

Bvals = [3.0800 0.0614 0.2317 0.8140 ...
    0.0130 2.3551 0.0486 0.0616 ...
    0.0482 0.5288 1.9321 0.0236 ...
    1.8563 0.0681 16.0480];


A = full(symtensor(Avals, 4,3));
B = full(symtensor(Bvals, 4,3));

%% Compute an eigenpair and test it
info = eig_geap(A, B, 'MaxIts', 100, 'Display',1);

%% 
% Test generalized eigenvector properties
x = info.x;
lambda = info.lambda;
lhs = ttsv(A,x,-1);
rhs = lambda * ttsv(B,x,-1);
fprintf('Generalized eigenpair identity norm: %g\n', norm(lhs-rhs));

%% Compute another eigenpair and test it
info = eig_geap(A, B, 'MaxIts', 100, 'Display',1);
%% 
% Test generalized eigenvector properties
x = info.x;
lambda = info.lambda;
lhs = ttsv(A,x,-1);
rhs = lambda * ttsv(B,x,-1);
fprintf('Generalized eigenpair identity norm: %g\n', norm(lhs-rhs));
