%% Symmetric CP Decomposition for Symmetric Tensors
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="cp.html">CP Decompositions</a> 
% &#62;&#62; <a href="cp_sym_doc.html">CP-SYM</a>
% </p>
% </html>
%
% The function |cp_sym| computes the symmetric CP decomposition of a
% symmetric tensor. 
% A *symmetric* tensor is invariant under any permutation
% of the indices. For a general dense |tensor|, we can verify its symmetry
% via the |issymmetric| function. An arbitrary dense tensor can be
% symmetrized by the |symmetrize| function. A symmetric tensor can also be
% stored as a |symtensor|. 
% The *symmetric CP decomposition* needs only a _single_ factor matrix
% (which is reused in every mode) and a weight/sign vector.  This can be
% stored as a |symktensor| object.
% The symmetric CP decompsition is described in the following reference:
%
% * T. G. Kolda, Numerical Optimization for Symmetric Tensor Decomposition,
% Mathematical Programming B, 151:225-248, 2015,
% <https://doi.org/10.1007/s10107-015-0895-0>  

%% Requirements
% Some of these codes requires an optimizaton solver to use. We recommend
% installing at least one of the following:
%
% * <https://github.com/sandialabs/poblano_toolbox *Poblano* Toolbox, Version 1.1>


%% Create a sample problem 
d = 3; % order
n = 10; % size
r = 2; % true rank

rng(5); % Set random number generator state for consistent results

info = create_problem('Size', n*ones(d,1), 'Num_Factors', r, ...
    'Symmetric', 1:d, 'Factor_Generator', @rand, 'Lambda_Generator', @rand, 'Noise', 0.1);

X = info.Data;
M_true = info.Soln; 
S_true = symktensor(M_true); % Convert from ktensor to symktensor

%%
% Check that the tensor is symmetric
issymmetric(X)

%%
% With even a small amount of noise, the gradient at the ideal solution can
% actually be large. This should be kept in mind when setting optimization
% termination conditions.

[f,g] = fg_explicit(S_true, X, norm(X)^2);
fprintf('Relative error || full(S_true) - X || / || X ||: %g\n', norm(full(M_true)-X)/norm(X));
fprintf('Function value at true solution: %g\n', f);
fprintf('Gradient norm at true solution: %g\n', norm(g));


%% Run CP-SYM using L-BFGS from Poblano Toolbox
% The default |cp_sym| uses a special objective function where each unique
% entry is only counted once in the sum-of-squared-error objective
% function. This really slows things down without much impact, so it's a
% good idea to set |'unique'| to |false|. Likewise, the |'l1param'|
% defaults to 10 but should be set to 0 for most problems. This is the
% recommended way to run the method:
optparams = lbfgs('defaults'); % Get the optimization parameters
optparams.RelFuncTol = 1e-10; % Tighten the stopping tolerance
optparams.StopTol = 1e-6; % Tighten the stopping tolerance
rng(5); % Set random number generator state for consistent results

[S,info] = cp_sym(X,r,'unique',false,'l1param',0,'alg_options',optparams);

fprintf('\n');
fprintf('Final function value: %.2g\n', fg_explicit(S, X, norm(X)^2));
fprintf('Stopping condition: %s\n', info.optout.ExitDescription);
fprintf('Check similarity score (1=perfect): %.2f\n', score(S,S_true));
fprintf('\n');

%% Run CP-SYM using FMINCON from Optimization Toolbox
% We can also run a version with nonnegativity constraints. In this case,
% we want to remove lambda from the optimization because otherwise there
% will be a scaling ambiguity. We need to use FMINCON because it accepts
% constraints.

rng(5); % Set random number generator state for consistent results

[S,info] = cp_sym(X,r,'unique',false,'l1param',0,'nonneg',true,'nolambda',true,'alg','fmincon');

fprintf('\n');
fprintf('Final function value: %.2g\n', fg_explicit(S, X, norm(X)^2));
fprintf('Check similarity score (1=perfect): %.2f\n', score(S,S_true));
%fprintf('Stopping condition: %s\n', info.optout.message);
fprintf('\n');
