function [xbest, fbest, info] = tt_opt_lbfgs(xinit, fgh, varargin)
%TT_OPT_LBFGS Wrapper for Poblano L-BFGS optimization method.
%
%   [X, F, INFO] = TT_OPT_LBFGS(X0, FGH, 'param', value, ...) is a wrapper
%   for the LBFGS method in Poblano. The wrapper just makes it a
%   bit easier to call from within Tensor Toolbox. Here X0 is in the
%   initial guess for the solution, FGH is a function handle to a function
%   that returns the function and gradient, and then there are optional
%   parameter-value pairs (see below).
%
%   The Poblano Toolbox is available on GITHUB:
%   https://github.com/sandialabs/poblano_toolbox/releases/tag/v1.2
%
%   For more optimzation algorithm choices and parameters, see
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','opt_options_doc.html')))">Tensor Toolbox Optimization Methods</a>, and
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','tt_opt_doc.html')))">Tensor Toolbox Optimization Methods for Developers</a>
%
%   See also TT_OPT_ADAM, TT_OPT_LBFGSB, TT_OPT_FMINUNC.
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
params.addParameter('printitn', 1); % DisplayIters
params.addParameter('m', 5); % M
params.addParameter('maxiters', 1000); % MaxIters
params.addParameter('subiters', 10) % MaxFuncEvals = maxiters * subiters
params.addParameter('ftol', 1e-10); % RelFuncTol
params.addParameter('gtol', 1e-5); % StopTol
params.addParameter('mdesc', 'Poblano L-BFGS Optimization');
params.addParameter('xdesc', []);
params.parse(varargin{:});

%% Options for printing
mdesc = params.Results.mdesc;
xdesc = params.Results.xdesc;

%% Setting optimization parameters
opts = params.Unmatched;
opts.M = params.Results.m;
opts.MaxIters = params.Results.maxiters;
opts.MaxFuncEvals = params.Results.maxiters * params.Results.subiters;
opts.RelFuncTol = params.Results.ftol;
opts.StopTol = params.Results.gtol;

printitn = params.Results.printitn;
if printitn == 0
    opts.Display = 'off';
else
    opts.DisplayIters = printitn;
end

opts.TraceFunc = true;
opts.TraceGradNorm = true;
opts.TraceRelFunc = true;

%% Welcome message
if printitn > 0
    fprintf('\n');
    fprintf('%s\n', mdesc);
    if isempty(xdesc)
        fprintf('Number of Variables: %d\n', n);
    else
        fprintf('%s\n', xdesc);
    end
    fprintf('Parameters: ');    
    fprintf('m=%d, ', opts.M);
    fprintf('ftol=%g, gtol=%g, ', opts.RelFuncTol, opts.StopTol);
    fprintf('maxiters = %d, maxsubiters=%d', opts.MaxIters, opts.MaxFuncEvals);
    fprintf('\n');
    fprintf('\n');
    fprintf('Begin Main Loop\n');
end
setuptime = toc(setupTimer);

%% Run optimization
tic; 
out = lbfgs(fgh, xinit, opts);
opttime = toc;

%% Save stuff
xbest = out.X;
fbest = out.F;
info.params = params.Results;
info.f_trace = out.TraceFunc;
info.gnorm_trace = out.TraceGradNorm;
info.time_trace = [];
info.setup_time = setuptime;
info.opt_time = opttime;
info.f_final = fbest;
info.iters = out.Iters;
info.subiters = out.FuncEvals;
info.exit_condition = out.ExitDescription;

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


info.out = out;

