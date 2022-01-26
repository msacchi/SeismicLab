%% Weighted Optimization for CP Tensor Decomposition with Incomplete Data
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="cp.html">CP Decompositions</a> 
% &#62;&#62; <a href="cp_wopt_doc.html">CP-WOPT</a>
% </p>
% </html>
%
% We explain how to use the CP Weighted Optimization (CP-WOPT) method
% implemented in |cp_wopt|. The method is described in the following article:
%
% * E. Acar, D. M. Dunlavy, T. G. Kolda and M. Mørup, 
%   Scalable Tensor Factorizations for Incomplete Data, 
%   Chemometrics and Intelligent Laboratory Systems, 106(1):41-56, 2011,
%   http://dx.doi.org/10.1016/j.chemolab.2010.08.004.

%% Third-party optimization software
% The |cp_wopt| method uses third-party optimization software to do the
% optimization. You can use either 
%
% * <https://github.com/stephenbeckr/L-BFGS-B-C *L-BFGS-B* by Stephen Becker> 
% (preferred), or
% * <https://software.sandia.gov/trac/poblano *POBLANO* Version 1.1 by
% Evrim Acar, Daniel Dunlavy, and Tamara Kolda>.
%
% The remainder of these instructions assume L-BFGS-B is being used. See
% <cp_wopt_poblano_doc.html here> for instructions on using |cp_wopt| with
% Poblano.

%% Important Information
% 
% It is critical to zero out the values in the missing entries of the data
% tensor. This can be done by calling |cp_wopt(X.*P,P,...)|. This is a
% frequent source of errors in using this method.

%% Create an example problem with missing data. 
% Here we have 25% missing data and 10% noise.   
R = 2;
info = create_problem('Size', [15 10 5], 'Num_Factors', R, ...
    'M', 0.25, 'Noise', 0.10);
X = info.Data;
P = info.Pattern;
M_true= info.Soln;

%% Create initial guess using 'nvecs'
M_init = create_guess('Data', X, 'Num_Factors', R, ...
    'Factor_Generator', 'nvecs');




%% Call the |cp_wopt| method
% Here is an example call to the cp_opt method. By default, each iteration
% prints the least squares fit function value (being minimized) and the
% norm of the gradient. 
[M,~,output] = cp_wopt(X, P, R, 'init', M_init);

%% Check the output
% It's important to check the output of the optimization method. In
% particular, it's worthwhile to check the exit message for any problems.
% The message |CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH| means that
% it has converged because the function value stopped improving.
exitmsg = output.ExitMsg


%% Evaluate the output
% We can "score" the similarity of the model computed by CP and compare
% that with the truth. The |score| function on ktensor's gives a score in
% [0,1]  with 1 indicating a perfect match. Because we have noise, we do
% not expect the fit to be perfect. See <matlab:doc('ktensor/score') doc
% score> for more details.
scr = score(M,M_true)

%% Create a SPARSE example problem with missing data. 
% Here we have 95% missing data and 10% noise.   
R = 2;
info = create_problem('Size', [150 100 50], 'Num_Factors', R, ...
    'M', 0.95, 'Sparse_M', true, 'Noise', 0.10);
X = info.Data;
P = info.Pattern;
M_true= info.Soln;

%% Create initial guess using 'nvecs'
M_init = create_guess('Data', X, 'Num_Factors', R, ...
    'Factor_Generator', 'nvecs');


%% Call the |cp_wopt| method
[M,~,output] = cp_wopt(X, P, R, 'init', M_init);

%% Check the output
exitmsg = output.ExitMsg


%% Evaluate the output
scr = score(M,M_true)

