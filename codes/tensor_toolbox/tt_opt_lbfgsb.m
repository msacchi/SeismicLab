function [xbest, fbest, info] = tt_opt_lbfgsb(xinit, fgh, varargin)
%TT_OPT_LBFGSB Wrapper for L-BFGS-B optimization code of Stephen Becker.
%
%   [X, F, INFO] = TT_OPT_LBFGSB(X0, FGH, 'param', value, ...) is a wrapper
%   for the LBFGSB wrapper by Stephen Becker. The wrapper just makes it a
%   bit easier to call from within Tensor Toolbox. Here X0 is in the
%   initial guess for the solution, FGH is a function handle to a function
%   that returns the function and gradient, and then there are options
%   parameter-value pairs as follows.
%
%   The lower and upper bounds are converted to vectors and passed as
%   arguments directly:
%      'lower' - Lower bounds, can be vector or scalar {-Inf}
%      'upper' - Upper bounds, can be vector or scalar {+Inf}
%
%   The following modify options for L-BFGS-B. The name of the option and
%   the default value that this wrapper uses are in curly braces. An
%   asterick means it is different than the L-BFGS-B default:
%      'maxiters' - Max outer iterations {maxIts, 1000*}
%      'printitn' - Printing frequency {printEvery, 1}
%      'm' - Limited-memory paramter {m, 5}
%      'subiters' - Max calls to FGH {maxTotalIts=maxiters*subiters, 10*}
%      'ftol' - Used as stopping condition {factr*eps, 1e-10}
%      'gtol' - Used as stopping condition {pgtol, 1e-5}
%
%   These options control the display:
%      'mdesc' - Method description, printed out {'L-BFGS-B Optimization'}
%      'xdesc' - Variable description {auto-generated}
%   
%   The L-BFGS-B code is available on GITHUB:
%   https://github.com/stephenbeckr/L-BFGS-B-C
%
%   For more optimzation algorithm choices and parameters, see
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','opt_options_doc.html')))">Tensor Toolbox Optimization Methods</a>, and
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','tt_opt_doc.html')))">Tensor Toolbox Optimization Methods for Developers</a>
%
%   See also TT_OPT_ADAM, TT_OPT_LBFGS, TT_OPT_FMINUNC.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%%
setupTimer = tic;

%% Setup
if ~iscolumn(xinit)
    error('Initial guess must be a column vector');
end

n = length(xinit);

%% Algorithm Parameters
params = inputParser;
params.KeepUnmatched = true;
params.addParameter('lower', -Inf); % 'l' 
params.addParameter('upper', +Inf); % 'u'
params.addParameter('maxiters', 1000); % maxIts
params.addParameter('printitn', 1); % printEvery
params.addParameter('m', 5); % m
params.addParameter('subiters', 10) % maxTotalIts = maxiters*subiters
params.addParameter('ftol', 1e-10);
params.addParameter('gtol', 1e-5); %pgtol
params.addParameter('mdesc', 'L-BFGS-B Optimization (https://github.com/stephenbeckr/L-BFGS-B-C)');
params.addParameter('xdesc', []);
params.parse(varargin{:});
mdesc = params.Results.mdesc;
xdesc = params.Results.xdesc;
printitn = params.Results.printitn;

%% Setup inputs for L-BFGS-B
opts = params.Unmatched; % Pass through any unmatched parameters
opts.x0 = xinit; 
opts.m = params.Results.m;
opts.factr = params.Results.ftol/eps;
opts.pgtol = params.Results.gtol;
opts.maxIts = params.Results.maxiters;
opts.maxTotalIts = params.Results.maxiters * params.Results.subiters;
opts.printEvery = printitn;
opts.errFcn = @(x) toc;

lb = params.Results.lower;
if isscalar(lb)
    lb = lb * ones(n,1);
elseif ~isequal(size(lb), size(xinit))
    error('Lower bound much be either a scalar or the same size as xinit');
end

ub = params.Results.upper;
if isscalar(ub)
    ub = ub * ones(n,1);
elseif ~isequal(size(ub), size(xinit))
    error('Upper bound much be either a scalar or the same size as xinit');
end

%% Welcome message
if printitn > 0
    fprintf('\n');
    fprintf('%s\n', mdesc);
    if isempty(xdesc)
        fprintf('Number of Variables: %d\n', n);
    else
        fprintf('%s\n', xdesc);
    end
    fprintf('Lower bound: ');
    if all(lb == lb(1))
        fprintf('%g', lb(1));
    else
        fprintf('user-specified vector');
    end
    fprintf(', ');
    fprintf('Upper bound: ');
    if all(ub == ub(1))
        fprintf('%g', ub(1));
    else
        fprintf('user-specified vector');
    end
    fprintf('\n');
    fprintf('Parameters: ');
    fprintf('m=%d, ', opts.m);
    fprintf('ftol=%g, ', opts.factr*eps);
    fprintf('gtol=%g, ', opts.pgtol);
    fprintf('maxiters = %d, ', opts.maxIts);
    fprintf('subiters = %d', opts.maxTotalIts/opts.maxIts);
    fprintf('\n');
    fprintf('\n');
    fprintf('Begin Main Loop\n');
end
setuptime = toc(setupTimer);

%% Run method
tic; % Do not rename this timer (used inside the lbfgsb method)
[xbest,fbest,optout] = lbfgsb(fgh, lb, ub, opts);
opttime = toc;

%% Save stuff
info.params = params.Results;
info.f_trace = optout.err(:,1);
info.gnorm_trace = optout.err(:,2);
info.time_trace = optout.err(:,3);
info.setup_time = setuptime;
info.opt_time = opttime;
info.f_final = fbest;
info.iters = optout.iterations;
info.subiters = optout.totalIterations;
info.exit_condition = optout.lbfgs_message1;

%% Goodbye message
if printitn > 0
    fprintf('End Main Loop\n');
    fprintf('\n');
    fprintf('Final f: %10.4e\n', info.f_final);
    fprintf('Setup time: %.2g seconds\n', info.setup_time);
    fprintf('Optimization time: %.2g seconds\n', info.opt_time);
    fprintf('Iterations: %d\n', info.iters);
    fprintf('Total iterations: %d\n', info.subiters);
    fprintf('Exit condition: %s\n', info.exit_condition);
end

