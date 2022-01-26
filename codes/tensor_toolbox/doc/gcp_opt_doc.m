%% Generalized CP (GCP) Tensor Decomposition
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="cp.html">CP Decompositions</a> 
% &#62;&#62; <a href="gcp_opt_doc.html">GCP-OPT</a>
% </p>
% </html>
%
% This document outlines usage and examples for the generalized CP (GCP)
% tensor decomposition implmented in |gcp_opt|. 
% GCP allows alternate objective functions besides sum of squared errors,
% which is the standard for CP.
% The code support both dense and sparse input tensors, but the sparse
% input tensors require randomized optimization methods. For some examples,
% see also <gcp_opt_amino_doc.html GCP-OPT and Amino Acids
% Dataset>.
%
% GCP is described in greater detail in the manuscripts: 
%
% * D. Hong, T. G. Kolda, J. A. Duersch, Generalized Canonical
%   Polyadic Tensor Decomposition, SIAM Review, 62:133-163, 2020,
%   <https://doi.org/10.1137/18M1203626>
% * T. G. Kolda, D. Hong, Stochastic Gradients for Large-Scale Tensor
%   Decomposition. SIAM J. Mathematics of Data Science, 2:1066-1095,
%   2020, <https://doi.org/10.1137/19m1266265>
%
%% Basic Usage
% The idea of GCP is to use alternative objective functions. As such, the
% most important thing to specify is the objective function. 
%
% The command |M = gcp_opt(X,R,'type',type)| 
% computes an estimate of the best rank-|R|
% generalized CP (GCP) decomposition of the tensor |X| for the specified
% generalized loss function specified by |type|. The input |X| can be a
% tensor or sparse tensor. The result |M| is a Kruskal tensor. 
% Some options for the objective function are:
%
% * |'binary'| - Bernoulli distribution for binary data
% * |'count'|  - Poisson distribution for count data (see also <cp_apr_doc.html |cp_apr|>)
% * |'normal'| - Gaussian distribution (see also <cp_als_doc.html |cp_als|> and <cp_opt_doc.html |cp_opt|>)
% * |'huber (0.25)'| - Similar to Gaussian but robust to outliers
% * |'rayleigh'| - Rayleigh distribution for nonnegative data
% * |'gamma'|  - Gamma distribution for nonnegative data 
%
% See <gcp_opt_fg_options_doc.html Function Types for GCP>
% for a complete list of options.

%% Manually specifying the loss function
% Rather than specifying a type, the user has the option to explicitly
% specify the objection function, gradient, and lower bounds using the
% following options:
%
% * |'func'| - Objective function handle, e.g., |@(x,m) (m-x).^2|
% * |'grad'| - Gradient function handle, e.g., |@(x,m) 2.*(m-x)|
% * |'lower'| - Lower bound, e.g., 0 or |-Inf|
%
% Note that the function must be able to work on vectors of x and m values.

%% Choice of Optimzation Method
% The default optimization method is L-BFGS-B (bound-constrained
% limited-memory BFGS). To use this, install the third-party software:
%
% * <https://github.com/stephenbeckr/L-BFGS-B-C *L-BFGS-B* by Stephen Becker> 
%
% The L-BFGS-B software can only be used for dense tensors.
% The other choice is to use a stochastic optimization method, either
% stochastic gradient descent (SGD) or ADAM. This can be used for
% dense or sparse tensors.
%
% The command |M = gcp_opt(X,R,...,'opt',opt)| specifies the optimization
% method where |opt| is one of the following strings:
%
% * |'lbfgsb'| - Bound-constrained limited-memory BFGS
% * |'sgd'| - Stochastic gradient descent (SGD)
% * |'adam'| - Momentum-based SGD method
%
% Each methods has parameters, which are described below.

%% Specifying Missing or Incomplete Data Using the Mask Option
% If some entries of the tensor are unknown, the method can mask off that
% data during the fitting process. To do so, specify a *mask* tensor |W|
% that is the same size as the data tensor |X|.  
% The mask tensor should be 1 if the entry in |X| is known and 0 otherwise.
% The call is |M = gcp_opt(X,R,'type',type','mask',W)|.

%% Other Options
% A few common options are as follows:
%
% * |'maxiters'| - Maximum number of outer iterations {1000}
% * |'init'| - Initialization for factor matrices {|'rand'|}
% * |'printitn'| - Print every n iterations; 0 for no printing {1}
% * |'state'| - Random state, to re-create the same outcome {[]}

%% Specifying L-BFGS-B Parameters
% In addition to the options above, there are two options used to modify
% the L-BFGS-B behavior. 
%
% * |'factr'| - Tolerance on the change on the objective value. Defaults to
% 1e7, which is multiplied by machine epsilon.
% * |'pgtol'| - Projected gradient tolerance, defaults to 1e-5.
%
% It can sometimes be useful to increase or decrease |pgtol| depending on
% the objective function and size of the tensor.


%% Specifying SGD and ADAM Parameters
% There are a number of parameters that can be adjusted for SGD and ADAM.
%
% *Stochastic Gradient.* There are three different sampling methods for
% computing the stochastic gradient:
%
% * _Uniform_ - Entries are selected uniformly at random. Default for dense
% tensors.
% * _Stratified_ - Zeros and nonzeros are sampled separately, which is
% recommended for sparse tensors. Default for sparse tensors.
% * _Semi-Stratified_ - Modification to stratified sampling that avoids
% rejection sampling for better efficiency at the cost of potentially
% higher variance.
%
% The options corresponding to these are as follows.
%
% * |'sampler'| - Type of sampling to use for stochastic gradient. Defaults
% to |'uniform'| for dense and |'stratified'| for sparse. The third options
% is |'semi-stratified'|.
% * |'gsamp'| - Number of samples for stochastic gradient. This should
% generally be O(sum(sz)*r). For the stratified or semi-stratified sampler,
% this can be two numbers. The first 
% is the number of nonzero samples and the second is the number of zero
% samples. If only one number is specified, then this is used as the number
% for both nonzeros and zeros, and the total number of samples is 2x what is
% specified.
%
% *Estimating the Function.* We also use sampling to estimate the function
% value. 
%
% * |'fsampler'| - This can be |'uniform'| (default for dense) or
% |'stratified'| (default for sparse) or a custom function handle.
% The custom function handleis primarily useful 
% in reusing the same sampled elements across different tests. For
% instance, we might create such a sampler by calling the hidden sampling
% function and saving its outputs:
%
%  [xsubs, xvals, wghts] = tt_sample_uniform(X, 10000);
%  fsampler = @() deal(xsubs, xvals, wghts);'
%
% * |'fsamp'| - Number of samples to estimate function.
% This should generally be somewhat large since we want this sample to
% generate a reliable estimate of the true function value. 
%
% The |'stratified'| sampler has an extra option: 
% * |'oversample'| - Factor to oversample when implicitly sampling zeros in
% the sparse case. Defaults to 1.1. Only adjust for very small tensors.
%
% There are some other options that are needed for SGD, the learning rate
% and a decrease schedule. Our schedule is very simple - we decrease the
% rate if there is no improvement in the approximate function value after
% an epoch. After a specified number of decreases (|'maxfails'|), we quit.
% 
% * |'rate'| - Initial learning rate. Defaults to 1e-3.
% * |'decay'| - How much to decrease the learning rate once progress
% stagnates, i.e., no decrease in objective function between epochs.
% Default to 0.1.
% * |'maxfails'| - How many times to decrease the learning rate. Can be
% zero. Defaults to 1.
% * |'epciters'| - Iterations per epoch. Defaults to 1000.
% * |'festtol'| - Quit if the function estimate goes below this level. Defaults to |-Inf|.
%
% There are some options that are specific to ADAM and generally needn't
% change:
%
% * |'beta1'| - Default to 0.9
% * |'beta2'| - Defaults to 0.999
% * |'epsilon'| - Defaults to 1e-8
% 

%% Example on Gaussian distributed 
% We set up the example with known low-rank structure. Here |nc| is the
% rank and |sz| is the size.
clear
rng('default')
nc = 2; 
sz = [50 60 70]; 
info = create_problem('Size',sz,'Num_Factors',nc);
X = info.Data;
M_true = info.Soln;
whos
%%
% Run GCP-OPT
tic, [M1,M0,out] = gcp_opt(X,nc,'type','normal','printitn',10); toc
fprintf('Final fit: %e (for comparison to f in CP-ALS)\n',1 - norm(X-full(M1))/norm(X));
fprintf('Score: %f\n',score(M1,M_true));

%%
% Compare to CP-ALS, which should usually be faster
tic, M2 = cp_als(X,nc,'init',tocell(M0),'printitn',1); toc
fprintf('Objective function: %e (for comparison to f(x) in GCP-OPT)\n', norm(X-full(M2))^2/prod(size(X)));
fprintf('Score: %f\n',score(M2,M_true));

%%
% Now let's try is with the ADAM functionality
tic, [M3,~,out] = gcp_opt(X,nc,'type','normal','opt','adam','init',M0,'printitn',1); toc
fprintf('Final fit: %e (for comparison to f in CP-ALS)\n',1 - norm(X-full(M1))/norm(X));
fprintf('Score: %f\n',score(M3,M_true));

%% Create an example Rayleigh tensor model and data instance.
% Consider a tensor that is Rayleigh-distribued. This means its entries are
% all nonnegative. First, we generate such a tensor with low-rank
% structure.
clear
rng(65)
nc = 3;
sz = [50 60 70];
nd = length(sz);

% Create factor matrices that correspond to smooth sinusidal factors
U=cell(1,nd);
for k=1:nd
    V = 1.1 + cos(bsxfun(@times, 2*pi/sz(k)*(0:sz(k)-1)', 1:nc));
    U{k} = V(:,randperm(nc));
end
M_true = normalize(ktensor(U));
X = tenfun(@raylrnd, full(M_true));
%%
% Visualize the true solution
viz(M_true, 'Figure', 1)

%%
% Run GCP-OPT
tic, [M1,~,out] = gcp_opt(X,nc,'type','rayleigh','printitn',10); toc
fprintf('Score: %f\n',score(M1,M_true));

%%
% Visualize the solution from GCP-OPT
viz(M1, 'Figure', 2)

%%
% Now let's try is with the scarce functionality - this leaves out all but
% 10% of the data!
tic, [M2,~,out] = gcp_opt(X,nc,'type','rayleigh','opt','adam'); toc
fprintf('Final fit: %e (for comparison to f in CP-ALS)\n',1 - norm(X-full(M1))/norm(X));
fprintf('Score: %f\n',score(M2,M_true));

%%
% Visualize the solution with scarce
viz(M2, 'Figure', 3)

%% Boolean tensor. 
% The model will predict the odds of observing a 1. Recall that the odds
% related to the probability as follows. If $p$ is the probability adn $r$
% is the odds, then $r = p / (1-p)$. Higher odds indicates a higher
% probability of observing a one. 
clear
rng(7639)
nc = 3; % Number of components
sz = [50 60 70]; % Tensor size
nd = length(sz); % Number of dimensions

%%
% We assume that the underlying model tensor has factor matrices with only
% a few "large" entries in each column. The small entries should correspond
% to a low but nonzero entry of observing a 1, while the largest entries,
% if multiplied together, should correspond to a very high likelihood of
% observing a 1.
probrange = [0.01 0.99]; % Absolute min and max of probabilities
oddsrange = probrange ./ (1 - probrange);
smallval = nthroot(min(oddsrange)/nc,nd);
largeval = nthroot(max(oddsrange)/nc,nd);

A = cell(nd,1);
for k = 1:nd
    A{k} = smallval * ones(sz(k), nc);
    nbig = 5;
    for j = 1:nc
        p = randperm(sz(k));
        A{k}(p(1:nbig),j) = largeval;
    end
end
M_true = ktensor(A);

%%
% Convert K-tensor to an observed tensor
% Get the model values, which correspond to odds of observing a 1
Mfull = full(M_true); 
% Convert odds to probabilities
Mprobs = Mfull ./ (1 + Mfull); 
% Flip a coin for each entry, with the probability of observing a one
% dictated by Mprobs 
Xfull = 1.0*(tensor(@rand, sz) < Mprobs); 
% Convert to sparse tensor, real-valued 0/1 tensor since it was constructed
% to be sparse
X = sptensor(Xfull); 
fprintf('Proportion of nonzeros in X is %.2f%%\n', 100*nnz(X) / prod(sz));

%%
% Just for fun, let's visualize the distribution of the probabilities in
% the model tensor. 
histogram(Mprobs(:))

%%
% Call GCP_OPT on the full tensor
[M1,~,out] = gcp_opt(Xfull, nc, 'type', 'binary','printitn',25);
fprintf('Final score: %f\n', score(M1,M_true));

%%
% GCP-OPT as sparse tensor

[M2,~,out] = gcp_opt(X, nc, 'type', 'binary');
fprintf('Final score: %f\n', score(M2,M_true));


%% Create and test a Poisson count tensor.
nc = 3;
sz = [80 90 100];
nd = length(sz);
paramRange = [0.5 60];
factorRange = paramRange.^(1/nd);
minFactorRatio = 95/100;
lambdaDamping = 0.8;
rng(21);
info = create_problem('Size', sz, ...
    'Num_Factors', nc, ...
    'Factor_Generator', @(m,n)factorRange(1)+(rand(m,n)>minFactorRatio)*(factorRange(2)-factorRange(1)), ...
    'Lambda_Generator', @(m,n)ones(m,1)*(lambdaDamping.^(0:n-1)'), ...
    'Sparse_Generation', 0.2);

M_true = normalize(arrange(info.Soln));
X = info.Data;
viz(M_true, 'Figure',3);

%% Loss function for Poisson negative log likelihood with identity link.

% Call GCP_OPT on sparse tensor
[M1,M0,out] = gcp_opt(X, nc, 'type', 'count','printitn',25);
fprintf('Final score: %f\n', score(M1,M_true));


