%% All-at-once optimization for CP tensor decomposition
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="cp.html">CP Decompositions</a> 
% &#62;&#62; <a href="cp_opt_doc.html">CP-OPT</a>
% </p>
% </html>
%
% We explain how to use |cp_opt| function which implements the *CP-OPT*
% method that fits the CP model using _direct_ or _all-at-once_
% optimization. This is in contrast to the |cp_als| function which
% implements the *CP-ALS* that fits the CP model using _alternating_ 
% optimization. The CP-OPT method is described in the
% following reference: 
%
% * E. Acar, D. M. Dunlavy and T. G. Kolda, A Scalable
% Optimization Approach for Fitting Canonical Tensor Decompositions,
% J. Chemometrics, 25(2):67-86, 2011,
% <http://doi.org/10.1002/cem.1335>


%% Third-party optimization software
% The |cp_opt| method uses third-party optimization software to do the
% optimization. You can use either 
%
% * <https://github.com/stephenbeckr/L-BFGS-B-C *L-BFGS-B* by Stephen Becker> 
% (preferred), or
% * <https://software.sandia.gov/trac/poblano *POBLANO* Version 1.1 by
% Evrim Acar, Daniel Dunlavy, and Tamara Kolda>.
%
% The remainder of these instructions assume L-BFGS-B is being used. See
% <cp_opt_poblano_doc.html here> for instructions on using |cp_opt| with
% Poblano.

%% Check that the software is installed. 
% Be sure that lbfgsb is in your path.
help lbfgsb

%% Create an example problem. 
% Create an example 50 x 40 x 30 tensor with rank 5 and add 10% noise.
R = 5;
info = create_problem('Size', [50 40 30], 'Num_Factors', R, 'Noise', 0.10);
X = info.Data;
M_true = info.Soln;

%% Create initial guess using 'nvecs'
M_init = create_guess('Data', X, 'Num_Factors', R, 'Factor_Generator', 'nvecs');

%% Call the |cp_opt| method
% Here is an example call to the cp_opt method. By default, each iteration
% prints the least squares fit function value (being minimized) and the
% norm of the gradient. 

[M,M0,output] = cp_opt(X, R, 'init', M_init);

%% Check the output
% It's important to check the output of the optimization method. In
% particular, it's worthwhile to check the exit message. 
% The message |CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH| means that
% it has converged because the function value stopped improving.
exitmsg = output.ExitMsg

%%
% The fit is the percentage of the data that is explained by the model.
% Because we have noise, we do not expect the fit to be perfect.
fit = output.Fit

%% Evaluate the output
% We can "score" the similarity of the model computed by CP and compare
% that with the truth. The |score| function on ktensor's gives a score in
% [0,1]  with 1 indicating a perfect match. Because we have noise, we do
% not expect the fit to be perfect. See <matlab:doc('ktensor/score') doc
% score> for more details.
scr = score(M,M_true)

%% Overfitting example
% Re-using the same example as before, consider the case where we don't
% know R in advance. We might guess too high. Here we show a case where we
% guess R+1 factors rather than R. 

% Generate initial guess of the corret size
M_plus_init = create_guess('Data', X, 'Num_Factors', R+1, ...
    'Factor_Generator', 'nvecs');

%%

% Run the algorithm
[M_plus,~,output] = cp_opt(X, R+1, 'init', M_plus_init);
exitmsg = output.ExitMsg
fit = output.Fit

%%

% Check the answer (1 is perfect)
scr = score(M_plus, M_true)

%% Nonnegative factorization
% We can employ lower bounds to get a nonnegative factorization.

%% Create an example problem. 
% Create an example 50 x 40 x 30 tensor with rank 5 and add 10% noise. We
% select nonnegative factor matrices and lambdas. The
% create_problem doesn't really know how to add noise without going
% negative, so we _hack_ it to make the observed tensor be nonzero.
R = 5;
info = create_problem('Size', [50 40 30], 'Num_Factors', R, 'Noise', 0.10,...
    'Factor_Generator', 'rand', 'Lambda_Generator', 'rand');
X = info.Data .* (info.Data > 0); % Force it to be nonnegative
M_true = info.Soln;

%% Generate initial guess of the corret size
M_init = create_guess('Data', X, 'Num_Factors', R, ...
    'Factor_Generator', 'rand');
%% Call the |cp_opt| method
% Here we specify a lower bound of zero with the last two arguments.
[M,M0,output] = cp_opt(X, R, 'init', M_init,'lower',0);

%% Check the output
exitmsg = output.ExitMsg

%%
% The fit is the percentage of the data that is explained by the model.
% Because we have noise, we do not expect the fit to be perfect.
fit = output.Fit

%% Evaluate the output
% We can "score" the similarity of the model computed by CP and compare
% that with the truth. The |score| function on ktensor's gives a score in
% [0,1]  with 1 indicating a perfect match. Because we have noise, we do
% not expect the fit to be perfect. See <matlab:doc('ktensor/score') doc
% score> for more details.
scr = score(M,M_true)
