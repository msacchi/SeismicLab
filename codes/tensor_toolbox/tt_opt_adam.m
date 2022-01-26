function [xbest, fbest, info] = tt_opt_adam(xinit, fh, gh, varargin)
%TT_OPT_ADAM ADAM stochastic optimization method.
%
%   [X,F,INFO] = TT_OPT_ADAM(X0, FH, GH, 'param', value, ...) applies the
%   Adam stochastic optimization method. This method depends on an
%   approximate/stochastic gradient, grouping iterations into epochs. Each
%   individual iteration may not yield improvement but should yield
%   improvement in aggreate. Therefore, iterations are grouped into
%   'epochs'. We track the process by evaluating the function value after
%   each epoch. The function values need not be exact. The user passes
%   function handles for evaluating the function and gradient. Here, X0 is
%   a column vector which is the intial guess, FH is a function handle such
%   that FH(X) evaluates the (approximate) function at a point X and
%   returns a scalar value, and GH is a function handle such that GH(X)
%   evaluates the stochastic gradient at a point X and returns a vector the
%   same size as X. It also takes a number of optional arguments to specify
%   algorithm parameters. These are specified below. The return values are
%   the vector X and function value F for the optimal function value. The
%   INFO return value includes extra information such as traces of the
%   function values over time, exit condition, algorithm  parameters, etc. 
%
%   The optional parameter-value pairs are as follows, with defaults in
%   curly braces. 
%
%   Key Algorithm Parameters
%      'lower' - Scalar lower bound on X {-Inf}
%      'subiters' - Number of iterations per epoch {100}
%      'maxiters' - Maximum number of epochs {100}
%      'rate' - Learning rate {1e-2}
%      'maxfails' - Max number of failed epochs allowed {1}
%
%   Minor Algorithm Paramters
%      'decay' - Decay of learning rate after failed epoch {0.1}
%      'backup' - Revert to end of previous epoch after failure {true}
%      'ftol' - Quit if function value goes below this value {-Inf}
%      'beta1', 'beta2', 'epsilon' - Adam parameters {0.9, 0.999, 1e-8}
%      'state' - State of random number generator {current state}
%
%   Output Control Parameters
%      'printitn' - Verbosity, 0 = none {1}
%      'mdesc' - Method description {'Adam Stochastic Optimization'}
%      'xdesc' - Description of input (e.g., size) {auto-generated}
%      'fdesc' - Description of (approximate) function computation {none}
%      'gdesc' - Description of stochastic gradient computation {none}
%      'fexact' - Boolean if function is computed exactly {true}
%   
%   References
%   * Original Method: D. P. Kingma and J. Ba, Adam: A Method for
%     Stochastic Optimization, December 2015, arXiv:1412.6980v9   
%   * Details on Stopping Conditions: T. G. Kolda and D. Hong,
%     Stochastic Gradients for Large-Scale Tensor Decomposition, SIAM
%     Journal on Mathematics of Data Science, accepted for publication,
%     arXiv:1906.01687   
%
%   For more optimzation algorithm choices and parameters, see
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','opt_options_doc.html')))">Tensor Toolbox Optimization Methods</a>, and
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','tt_opt_doc.html')))">Tensor Toolbox Optimization Methods for Developers</a>
%
%   See also TT_OPT_LBFGSB, TT_OPT_LBFGS, TT_OPT_FMINUNC.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

%%
setupTimer = tic;

%% Random set-up
defaultStream = RandStream.getGlobalStream;

%% Algorithm Parameters
params = inputParser;
params.addParameter('lower', [], @isscalar);
params.addParameter('maxiters',100);
params.addParameter('printitn',1);
params.addParameter('state', defaultStream.State);
params.addParameter('rate',1e-2);
params.addParameter('decay',0.1);
params.addParameter('maxfails',1);
params.addParameter('subiters', 100);
params.addParameter('backup', true);
params.addParameter('ftol', -Inf);
params.addParameter('beta1', 0.9);
params.addParameter('beta2', 0.999);
params.addParameter('epsilon', 1e-8);
params.addParameter('mdesc', 'Adam Stochastic Optimization');
params.addParameter('xdesc', []);
params.addParameter('fdesc', []);
params.addParameter('gdesc', []);
params.addParameter('fexact', true);
params.parse(varargin{:});

%% Initialize random number generator with specified state
defaultStream.State = params.Results.state;

%% Setup
lb = params.Results.lower;
maxiters = params.Results.maxiters;
printitn = params.Results.printitn;
rate = params.Results.rate;
decay = params.Results.decay;
maxfails = params.Results.maxfails;
subiters = params.Results.subiters;
backup = params.Results.backup;
festtol = params.Results.ftol;
beta1 = params.Results.beta1;
beta2 = params.Results.beta2;
epsilon = params.Results.epsilon;
mdesc = params.Results.mdesc;
xdesc = params.Results.xdesc;
fdesc = params.Results.fdesc;
gdesc = params.Results.gdesc;
fexact = params.Results.fexact;

%% Error checking
if ~iscolumn(xinit)
    error('Initial guess must be a column vector');
end

if fexact
    fstr = 'f';
else
    fstr = 'f~';
end

%% Lower bounds
if isempty(lb)
    lb = -Inf;
end

%% Problem size
n = length(xinit);

%% Welcome message
if printitn > 0
    fprintf('\n');
    fprintf('%s\n', mdesc); % Method Description
    if isempty(xdesc)
        fprintf('Number of Variables: %d\n', n);
    else
        fprintf('%s\n', xdesc);
    end
    if lb ~= -Inf, fprintf('Variable lower bound: %g\n', lb); end
    if ~isempty(fdesc), fprintf('%s\n', fdesc); end
    if ~isempty(gdesc), fprintf('%s\n', gdesc); end
    fprintf('Max iterations (epochs): %d\n',maxiters);
    fprintf('Iterations per epoch: %d\n', subiters);
    fprintf('Learning rate / decay / maxfails: %g %g %g\n', rate, decay, maxfails);
    fprintf('Backup to best known solution after epoch fail? ');
    if backup, fprintf('true'); else, fprintf('false'); end
    fprintf('\n');
    fprintf('beta1 / beta2 / epsilon: %g %g %g\n', beta1, beta2, epsilon);
    fprintf('\n');
    fprintf('Begin Main Loop\n');
end

%% Tracing setup
fest_trace = zeros(maxiters+1,1);
step_trace = zeros(maxiters+1,1);
time_trace = zeros(maxiters+1,1);

%% Initialization

x = xinit;          % Initial Guess
fbest = +Inf;       % Best previous function value
step = 0;           % Last step
nfails = 0;         % Counter # times fails to decrease
titers = 0;         % Total iterations (used by ADAM to determine step)
m = zeros(n,1);     % Adam variable
v = zeros(n,1);     % Adam variable
setuptime = toc(setupTimer);

%% Main loop
% First iteration as recorded in the trace doesn't really count. This is
% stripped out when we save the final information.
mainTimer = tic;
for nepoch = 1:maxiters+1
    
    % Compute (possibly approximate) objective function value
    fest = fh(x);
    
    % Check for improvement in this epoch
    good_epoch = fest < fbest;
    
    % Increment number of failed epochs
    if ~good_epoch
        nfails = nfails + 1;
    end

    % Check various termination conditions
    if nfails > maxfails
        exit_now = true;
        exit_condition = sprintf('Number of failed epochs > %d', maxfails);
    elseif fest < festtol
        exit_now = true;
        exit_condition = sprintf('Function value < %g', festtol);
    elseif nepoch == maxiters+1
        exit_now = true;
        exit_condition = sprintf('Number of epochs > %d', maxiters);
    else
        exit_now = false;
    end
    
    % Reporting
    if (printitn > 0) && (mod(nepoch-1,printitn) == 0 || ~good_epoch || exit_now)
        fprintf('Epoch %2d: %s = %e', nepoch-1, fstr, fest);
        if nepoch > 1
            fprintf(', step = %g', step);
        end
        if ~good_epoch
            fprintf(', nfails = %d', nfails);
            if backup
                fprintf(' (resetting to solution from last epoch)');
            end
        end
        fprintf('\n');
        
    end
    
    if good_epoch % Save best solution so far        
        xbest = x;
        fbest = fest;        
        msave = m;
        vsave = v;
    elseif backup % Back up to best solution so far!
        x = xbest;
        m = msave;
        v = vsave;
        titers = titers - subiters;
    end
    
    % Save trace
    fest_trace(nepoch) = fest;
    step_trace(nepoch) = step;
    time_trace(nepoch) = toc(mainTimer);
    
    % Check if it's time to quit
    if exit_now
        break;
    end
    
    % Main loop inner iteration
    step = decay^nfails * rate;     
    
    for iter = 1:subiters
        
        % Tracking total iterations
        titers = titers + 1;
        
        % Compute estimated gradient
        gest = gh(x);
        
        % Check for inf
        isinfgrad = any(isinf(gest));
        if any(isinfgrad)
            error('Infinite gradient encountered! (epoch = %g, iter = %g)', nepoch, iter);
        end
        
        % Take a step
        m = beta1*m + (1-beta1)*gest;
        v = beta2*v + (1-beta2)*gest.^2;
        mhat = m/(1-beta1^titers);
        vhat = v/(1-beta2^titers);
        x = max(lb,x-step*mhat./(sqrt(vhat)+epsilon));
    end       
    
end

%% Finished, clean up
info.params = params.Results;
info.f_trace = fest_trace(2:nepoch);
info.step_trace = step_trace(2:nepoch);
info.time_trace = time_trace(2:nepoch) - time_trace(1);
info.f_init = fest_trace(1);
info.setup_time = setuptime + time_trace(1);
info.opt_time = info.time_trace(end);
info.f_final = fbest;
info.epochs = nepoch-1;
info.iters = info.epochs * subiters;
info.exit_condition = exit_condition;

%% Goodbye message
if printitn > 0
    fprintf('End Main Loop\n');
    fprintf('\n');
    fprintf('Final %s: %10.4e\n', fstr, info.f_final);
    fprintf('Setup time: %.2f seconds\n', info.setup_time);
    fprintf('Optimization time: %.2f seconds\n', info.opt_time);
    fprintf('Number of epochs: %d\n', info.epochs);
    fprintf('Total iterations: %d\n', info.iters);
    fprintf('Exit condition: %s\n', info.exit_condition);
end

