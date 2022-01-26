%% Developer Information for Optimization Methods in Tensor Toolbox
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="opt_options_doc.html">Optimization Methods</a> 
% &#62;&#62; <a href="tt_opt_doc.html">Developer Information</a> 
% </p>
% </html>
%
% This page contains information for developers who want to use the
% optimization method interfaces contained in Tensor Toolbox. These
% interfaces are meant to standardize across packages with very different
% ways of calling the functions, different parameter names, etc.
% For optimzation algorithm installation instructions and parameter
% options, see <opt_options_doc.html 
% Documentation for Tensor Toolbox Optimization Methods>.
%
% The above-referenced document is for general users. The information below
% is for developers of Tensor Toolbox methods that want to use these
% standardized optimization algorithms.

%% Test optimization problem 
% We use a simple but notoriously difficult test function to demonstrate
% the optimization methods, a variation on Rosenbrock:
%
% $$f(x) = (x_1 - 1)^2 + 4 \sum_{i=2}^n (x_i - x_{i-1}^2)^2$$
%
% This function is implemented in |optfunc| which takes |true| as a second
% argument to return only the gradient. The number of parameters (n) is
% automatically determined by the method. The optimum is the all-ones 
% vector. The code for the function is as follows:
%
% <include>optfunc.m</include>
%

%% 
% We create function handles to use the optimization solvers.

% Create a function handle to evaluate both f and g
fgh = @optfunc;
% Create a function handle to evaluate just f
fh = @optfunc;
% Create a function handle to evaluate just g
gh = @(x) optfunc(x,true);

%%
% Set up initial guess and optimum
n = 10;
xinit = [2; ones(n-1,1)];
xopt = ones(n,1);

% Initial function value
finit = fh(xinit);
fprintf('Function value at xinit = %g\n', finit);

% Optimal function value
fopt = fh(xopt);
fprintf('Function value at xopt = %g\n', fopt);

%% Sandia Poblano Toolbox L-BFGS 
% <https://software.sandia.gov/trac/poblano *POBLANO* Version 1.1 by
% Evrim Acar, Daniel Dunlavy, and Tamara Kolda>
[x,f,info] = tt_opt_lbfgs(xinit, fgh);

%%
% Check how close to optimal
fprintf('||x-xopt||/||xopt|| = %e\n', norm(x-xopt)/norm(xopt));
%%
% Tighten it up a bit by changing the tolerances. ALso note that we changed
% the printing frequency to once every *5* iterations rather than every
% iteration as above.
[x,f,info] = tt_opt_lbfgs(xinit, fgh, 'printitn', 5, 'gtol', 1e-20, 'ftol', 1e-25);

%%
% With tighter tolerances, we get closer to optimal
fprintf('||x-xopt||/||xopt|| = %e\n', norm(x-xopt)/norm(xopt));

%% Stephen Becker L-BFGS-B 
% <https://github.com/stephenbeckr/L-BFGS-B-C *L-BFGS-B* by Stephen Becker>.  
[x,f,info] = tt_opt_lbfgsb(xinit, fgh);

%%
% Check how close to optimal
fprintf('||x-xopt||/||xopt|| = %e\n', norm(x-xopt)/norm(xopt));


%% 
% Tighten it up a bit as before
[x,f,info] = tt_opt_lbfgsb(xinit, fgh, 'printitn', 5, 'gtol', 1e-12, 'ftol', 1e-24);

%%
% With tigher tolerances, we again get closer to optimal
fprintf('||x-xopt||/||xopt|| = %e\n', norm(x-xopt)/norm(xopt));

%% MATLAB Optimization Toolbox |fminunc| (Quasi-Newton Method)
% Distributed with the MATLAB Optimization Toolbox.
[x,f,info] = tt_opt_fminunc(xinit, fgh);

%%
% Check how close to optimal
fprintf('||x-xopt||/||xopt|| = %e\n', norm(x-xopt)/norm(xopt));

%% 
% Tighten it up a bit as before. We don't have control over the print
% frequency with this one, nor the change in function tolerance. But we can
% pass through the `StepTolerance`.
[x,f,info] = tt_opt_fminunc(xinit, fgh, 'gtol', 1e-12, 'StepTolerance',1e-14);

%%
% With tigher tolerances, we again get closer to optimal
fprintf('||x-xopt||/||xopt|| = %e\n', norm(x-xopt)/norm(xopt));

%% Adam (internal to Tensor Toolbox)
% Adam is a stochastic optimization method, so it's not especially fair to
% compare it to the other methods. But we do it anyways.
rng('default')
[x,f,info] = tt_opt_adam(xinit, fh, gh,'subiters',500,'maxfails',2,'maxiters',20000,'printitn',500,'fdesc','Function: exact','gdesc','Gradient: exact');
%%
% Check how close to optimal
fprintf('||x-xopt||/||xopt|| = %e\n', norm(x-xopt)/norm(xopt));

%%
% Let's tighten the convergence tolerance.
rng(1)
[x,f,info] = tt_opt_adam(xinit, fh, gh,'subiters',500,'maxfails',6,'maxiters',20000,'printitn',500,'fdesc','Function: exact','gdesc','Gradient: exact');
%%
% Check how close to optimal
fprintf('||x-xopt||/||xopt|| = %e\n', norm(x-xopt)/norm(xopt));

