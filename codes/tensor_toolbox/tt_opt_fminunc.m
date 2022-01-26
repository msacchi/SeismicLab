function [xbest, fbest, info] = tt_opt_fminunc(xinit, fgh, varargin)
%TT_OPT_FMINUNC Wrapper for FMINUNC in MATLAB Optimization Toolbox.
%
%   [X, F, INFO] = TT_OPT_FMINUNC(X0, FGH, 'param', value, ...) is a
%   wrapper for the Quasi-Newton method in the MATLAB Optimization Toolbox.
%   The wrapper just makes it a bit easier to call from within Tensor
%   Toolbox. Here X0 is in the initial guess for the solution, FGH is a
%   function handle to a function that returns the function and gradient,
%   and then there are optional  parameter-value pairs (see below).
%
%   For more optimzation algorithm choices and parameters, see
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','opt_options_doc.html')))">Tensor Toolbox Optimization Methods</a>, and
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','tt_opt_doc.html')))">Tensor Toolbox Optimization Methods for Developers</a>
%
%   See also TT_OPT_ADAM, TT_OPT_LBFGSB, TT_OPT_LBFGS.
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
params.addParameter('maxiters', 1000); % maxIts
params.addParameter('printitn', 1); % printEvery
params.addParameter('subiters', 10) % maxTotalIts = maxiters*subiters
params.addParameter('gtol', 1e-5); %pgtol
params.addParameter('mdesc', 'Quasi-Newton Optimization (via Optimization Toolbox)');
params.addParameter('xdesc', []);
params.parse(varargin{:});
mdesc = params.Results.mdesc;
xdesc = params.Results.xdesc;
printitn = params.Results.printitn;

%% Setup options for quasi-Newton method
opts = optimoptions('fminunc');
opts.SpecifyObjectiveGradient = true;
opts.MaxFunctionEvaluations = params.Results.maxiters * params.Results.subiters;
opts.MaxIterations = params.Results.maxiters;
opts.OptimalityTolerance = params.Results.gtol;
if printitn == 0
    opts.Display = 'off';
else
    opts.Display = 'iter';
end

moreparams = fieldnames(params.Unmatched);
if ~isempty(moreparams)
    for i = 1:length(moreparams)
       opts.(moreparams{i}) = params.Unmatched.(moreparams{i}); 
    end   
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
    fprintf('Parameters: ');
    fprintf('gtol=%g, ', opts.OptimalityTolerance);
    fprintf('maxiters = %d, maxsubiters=%d', opts.MaxIterations, opts.MaxFunctionEvaluations);
    fprintf('\n');
    fprintf('\n');
    fprintf('Begin Main Loop\n');
end
setuptime = toc(setupTimer);

%% Run method
tic;
[xbest, fbest, exitflag, output] = fminunc(fgh, xinit, opts);
opttime = toc;

%% Save stuff
info.params = params.Results;
info.exitflag = exitflag;
info.setup_time = setuptime;
info.opt_time = opttime;
info.f_final = fbest;
info.iters = output.iterations;
info.subiters = output.funcCount;
info.opts = opts;
switch(exitflag)
    case 1
        info.exit_condition = 'Gradient norm smaller than tolerance';
    case 2
        info.exit_condition = 'Change in x smaller than tolerance';
    case 3
        info.exit_condition = 'Change in function smaller than tolerance';
    case 5
        info.exit_condition = 'Line search failure';
    case 0
        info.exit_condition = 'Iterations exceeded max allowed';
    otherwise
        info.exit_condition = 'Unexpected failure';
end

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


info.opt_time = opttime;
info.f_final = fbest;
info.output = output;
