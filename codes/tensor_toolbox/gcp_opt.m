function [M, M0, info] = gcp_opt(X, nc, varargin)
%GCP_OPT Fits Generalized CP decomposition with user-specified function.
%
%   M = GCP_OPT(X,R,'type',TYPE) computes an estimate of the best rank-R
%   generalized CP (GCP) decomposition of the tensor X for the specified
%   generalized loss function specified by TYPE. The input X can be a
%   tensor or sptensor. The result M is a ktensor. A full set of types can
%   be found in GCP_FG_SETUP, but popular options include    
%
%      'binary' - Bernoulli distribution for binary data
%      'count'  - Poisson distribution for count data (see also CP_APR)
%      'normal' - Gaussian distribution (see also CP_ALS and CP_OPT)
%      'huber (0.25)' - Similar to Gaussian but robust to outliers
%      'rayleigh' - Rayleigh distribution for nonnegative data
%      'gamma'  - Gamma distribution for nonnegative data 
%
%   M = GCP_OPT(X,R,'func',FH,'grad',GH,'lower',LB) passes a user-specified
%   choice for the elementwise loss function, corresponding gradient, and
%   lower bound on the factor matrix entries. This is an alternative to
%   specifying the 'type' as shown above.
%
%   M = GCP_OPT(X,R,...,'opt',OPT) specifies the optimization method:
%      'lbfgsb' - Bound-constrained limited-memory BFGS
%      'sgd' - Stochastic gradient descent (SGD)
%      'adam' - Momentum-based SGD method
%   If X is dense, any of the three options can be used, and 'lbfgsb' is
%   the default. The X is sparse, only 'sgd' and 'adam' (default) are
%   options. Each method has specific parameters, see the documentation for
%   details.
%
%   M = GCP_OPT(X,R,...,'mask',W) specifies a mask W that is 0 for missing
%   entries and 1 everywhere else. This can only be used in the case that X
%   is dense and 'opt' is the default ('lbfgsb'). The missing entries are
%   ignored in the fitting of the model.
%
%   M = GCP_OPT(X,R,'param',value,...) specifies additional optional
%   parameters and values as follows, with the defaults in curly braces:
%
%      'maxiters' - Maximum number of outer iterations {1000}
%      'init' - Initialization for factor matrices {'rand'}
%      'printitn' - Print every n iterations; 0 for no printing {1}
%      'state' - Random state, to re-create the same outcome {[]}
%
%   [M,M0,out] = GCP_OPT(...) also returns the initial guess (M0) and a
%   structure with additional information. To reproduce the
%   run exactly, use M_alt = gcp_opt(X,R,out.params.Results).
%
%   REFERENCES: 
%   * D. Hong, T. G. Kolda, J. A. Duersch, Generalized Canonical
%     Polyadic Tensor Decomposition, SIAM Review, 62:133-163, 2020,
%     https://doi.org/10.1137/18M1203626     
%   * T. G. Kolda, D. Hong, Stochastic Gradients for Large-Scale Tensor
%     Decomposition. SIAM J. Mathematics of Data Science, 2:1066-1095,
%     2020, https://doi.org/10.1137/19m1266265
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','gcp_opt_doc.html')))">Documentation page for GCP_OPT</a>
%
%   See also CP_OPT, CP_APR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

% Created by Tamara G. Kolda, Fall 2018. Includes work with
% collaborators David Hong and Jed Duersch. 

%% Timers
setupStartA = tic;

%% Iniital setup
nd  = ndims(X);
sz = size(X);
tsz = prod(sz);

%% Random set-up
defaultStream = RandStream.getGlobalStream;

%% Set algorithm parameters from input or by using defaults
params = inputParser;

params.addParameter('type', [], @ischar);
params.addParameter('func', [], @(x) isa(x,'function_handle'));
params.addParameter('grad', [], @(x) isa(x,'function_handle'));
params.addParameter('lower', [], @isscalar);

params.addParameter('opt', [], @ischar);
params.addParameter('mask', [], @(x) isa(x,'tensor'));

params.addParameter('maxiters', 1000, @isscalar);
params.addParameter('init','rand');
params.addParameter('printitn',1, @isscalar);
params.addParameter('state', defaultStream.State);

params.addParameter('factr', 1e7, @isscalar);
params.addParameter('pgtol', 1e-4 * tsz, @isscalar);

params.addParameter('fsamp',[]);
params.addParameter('gsamp',[]);
params.addParameter('oversample',1.1, @isscalar);
params.addParameter('sampler', []);
params.addParameter('fsampler', []);

params.addParameter('rate', 1e-3, @isscalar);
params.addParameter('decay', 0.1, @isscalar);
params.addParameter('maxfails', 1, @isscalar);
params.addParameter('epciters', 1000, @isscalar);
params.addParameter('festtol', -Inf, @isscalar);

params.addParameter('beta1', 0.9, @isscalar);
params.addParameter('beta2', 0.999, @isscalar);
params.addParameter('epsilon', 1e-8, @isscalar);


params.parse(varargin{:});

% Save info
info.params = params.Results;

%% Initialize random number generator with specified state
defaultStream.State = params.Results.state;

%% Extract remaining parameters
type = params.Results.type;
fh = params.Results.func;
gh = params.Results.grad;
lb = params.Results.lower;
opt = params.Results.opt;
W = params.Results.mask;
maxiters = params.Results.maxiters;
init = params.Results.init;
printitn = params.Results.printitn;
factr = params.Results.factr;
pgtol = params.Results.pgtol;
fsamp = params.Results.fsamp;
gsamp = params.Results.gsamp;
oversample = params.Results.oversample;
gsampler_type = params.Results.sampler;
fsampler_type = params.Results.fsampler;
rate = params.Results.rate;
decay = params.Results.decay;
maxfails = params.Results.maxfails;
epciters = params.Results.epciters;
festtol = params.Results.festtol;
beta1 = params.Results.beta1;
beta2 = params.Results.beta2;
epsilon = params.Results.epsilon;

%% More setup
vecsz = sum(sz)*nc;
isdense = isa(X,'tensor');
issparse = isa(X,'sptensor');

if isdense  
    if isempty(W)
        nmissing = 0;
        nnonzeros = nnz(X);
        nzeros = tsz - nnonzeros;
    else
        X = X .* W;
        nmissing = tsz - nnz(W); 
        nnonzeros = nnz(X);
        nzeros = nnz(W) - nnz(X);
    end    
elseif issparse    
    if ~isempty(W)
        error('Cannot specify missing entries for sparse tensors');
    end
    nmissing = 0;
    nnonzeros = nnz(X);
    nzeros = tsz - nnonzeros;    
else    
    error('Input tensor must be tensor or sptensor');    
end

% Save info
info.tsz = tsz;
info.nmissing = nmissing;
info.nnonzeros = nnonzeros;
info.nzeros = nzeros;

%% Set up function/gradient and bounds

if ~isempty(fh) && ~isempty(gh) 
    if isempty(fh) || isempty(gh) 
        error('Must specify ''func'' and ''grad'' if either one is specified');
    end
    if isempty(lb)
        lb = -Inf;
    end
    type = 'user-specified';
else     
    [fh,gh,lb_] = tt_gcp_fg_setup(type,X);
    if isempty(lb)
        lb = lb_;
    end
end

% Save info
info.type = type;
info.fh = fh;
info.gh = gh;
info.lb = lb;

%% Create initial guess, denoted M0
if iscell(init)
    Uinit = init;
    M0 = ktensor(Uinit);
    inittype = 'cell';
elseif isa(init, 'ktensor')
    M0 = init;
    inittype = 'ktensor';
elseif strcmp(init,'rand')
    Uinit = cell(nd,1);
    for k = 1:nd
        Uinit{k} = rand(sz(k),nc);
    end
    M0 = ktensor(Uinit);
    M0 = M0 * (norm(X)/norm(M0)); % normalize
    inittype = 'rand';
end

% We assume throughout that the lambda weights are all ones. Make sure that
% the initial guess satisfies this property.
M0 = normalize(M0,0); 

% Save info
info.inittype = inittype;

%% Setup the optimization

if isempty(opt)
    if isdense
        opt = 'lbfgsb';
    elseif issparse
        opt = 'adam';
    end
end


if ~ismember(opt,{'lbfgsb','sgd','adam','adagrad'})
    error('Invalid choice for ''opt''');
end

use_lbfgsb = strcmpi(opt,'lbfgsb');
use_adam = strcmpi(opt,'adam');
use_sgd = strcmpi(opt,'sgd');
use_adagrad = strcmpi(opt,'adagrad');

use_stoc = use_adam || use_sgd || use_adagrad;

if issparse && ~use_stoc
    error('Must set ''opt'' to ''sgd'' or ''adam'' or ''adagrad'' for sparse tensor');
end

% Save info
info.opt = opt;

%% Set up for Stochastic Optimization
if use_stoc
    
    if ~isempty(W)
        error('Have not yet implemented stochastic optimization for the case of missing data');
    end
    
    crng = []; % Default value
    xnzidx = []; % Only create the sorted indices if needed
    
    
    % Set up fsampler
    if isempty(fsampler_type)
        if issparse
            fsampler_type = 'stratified';
        else
            fsampler_type = 'uniform';
        end
    end
        
    if isa(fsampler_type,'function_handle')
        
        fsampler = fsampler_type;
        fsampler_str = 'user-specified';
        
    elseif strcmp(fsampler_type, 'stratified')
        
        if isempty(fsamp)
            ftmp = max(ceil(nnonzeros/100), 10^5);
            fsamp(1) = min(ftmp, nnonzeros);
            fsamp(2) = min([ftmp, nnonzeros, nzeros]);
        elseif length(fsamp) == 1
            tmp = fsamp;
            fsamp(1) = tmp;
            fsamp(2) = tmp;
        end
        
        % Create and sort linear indices of X nonzeros for the sampler
        if isempty(xnzidx)
            xnzidx = tt_sub2ind64(sz,X.subs);
            xnzidx = sort(xnzidx);
        end
        
        fsampler = @() tt_sample_stratified(X, xnzidx, fsamp(1), fsamp(2), oversample);
        fsampler_str =  sprintf('stratified with %d nonzero and %d zero samples', fsamp);
        
        
    elseif strcmp(fsampler_type, 'uniform')
        
        if isempty(fsamp)
            fsamp = min( max(ceil(tsz/10), 10^6), tsz );
        end
        
        fsampler = @() tt_sample_uniform(X,fsamp);
        fsampler_str = sprintf('uniform with %d samples', fsamp);
        
    else
        
        error('Invalid choice for ''fsampler''');
        
    end
    
    % Set up gsampler
    if isempty(gsampler_type)
        if issparse 
            gsampler_type = 'stratified';
        else
            gsampler_type = 'uniform';
        end
    end
    
    if strcmp(gsampler_type, 'semi-stratified') || strcmp(gsampler_type, 'stratified')
        
        if isempty(gsamp)
            gtmp = max(1000, ceil(3 * nnonzeros / maxiters));
            gsamp(1) = min(gtmp, nnonzeros);
            gsamp(2) = min([gtmp, nnonzeros, nzeros]);
        end
        
        if length(gsamp) == 1
            tmp = gsamp;
            gsamp(1) = tmp;
            gsamp(2) = tmp;
        end
               
        if strcmp(gsampler_type, 'semi-stratified')
            
            gsampler = @() tt_sample_semistrat(X, gsamp(1), gsamp(2));
            gsampler_str = sprintf('semi-stratified with %d nonzero and %d zero samples', gsamp);
            crng = 1:gsamp(1);
        
        else      
            
            % Create and sort linear indices of X nonzeros for the sampler
            if isempty(xnzidx)
                xnzidx = tt_sub2ind64(sz,X.subs);
                xnzidx = sort(xnzidx);
            end

            gsampler = @() tt_sample_stratified(X, xnzidx, gsamp(1), gsamp(2), oversample);
            gsampler_str = sprintf('stratified with %d nonzero and %d zero samples', gsamp);
        end
        
    elseif strcmp(gsampler_type, 'uniform')     
        
        if isempty(gsamp)
            gsamp = min( max( 1000, ceil(10 * tsz / maxiters) ), tsz);
        end
        
        if issparse
            
            exp_nonzeros = gsamp * nnonzeros / tsz;
            exp_zeros = gsamp * nzeros / tsz;
            
            % Create and sort linear indices of X nonzeros for the sampler
            if isempty(xnzidx)
                xnzidx = tt_sub2ind64(sz,X.subs);
                xnzidx = sort(xnzidx);
            end            
            
            gsampler = @() tt_sample_stratified(X, xnzidx, random('Poisson', exp_nonzeros), random('Poisson', exp_zeros), oversample);
            gsampler_str = sprintf('pseudo-uniform with %d samples', gsamp);
            
            
        else
            
            gsampler = @() tt_sample_uniform(X,gsamp);
            gsampler_str = sprintf('uniform with %d samples', gsamp);
            
        end
    else
        
        error('Invalid sampler: %s', gsampler_type);
    
    end
    

    % Save info
    info.fsampler = fsampler_str;
    info.gsampler = gsampler_str;
    info.fsamp = fsamp;
    info.gsamp = gsamp;

end

setupTimeA = toc(setupStartA);

%% Welcome message
if printitn > 0
    fprintf('\n');
    fprintf('GCP-OPT-%s (Generalized CP Tensor Decomposition)\n',upper(opt));
    fprintf('\n');
    fprintf('Tensor size: %s (%d total entries)\n', tt_size2str(size(X)), tsz);
    if nmissing > 0
        fprintf('Missing entries: %d (%.2g%%)\n', nmissing, 100*nmissing / tsz);
    end
    if issparse
        fprintf('Sparse tensor: %d (%.2g%%) Nonzeros and %d (%.2f%%) Zeros\n', nnonzeros, 100*nnonzeros/tsz, nzeros, 100*nzeros/tsz);
    end
    fprintf('Generalized function Type: %s\n', type);
    fprintf('Objective function: %s\n', func2str(fh));
    fprintf('Gradient function: %s\n', func2str(gh));
    fprintf('Lower bound of factor matrices: %g\n', lb);
    fprintf('Optimization method: %s\n', opt);
    if use_stoc 
        fprintf('Max iterations (epochs): %d\n',maxiters);
        fprintf('Iterations per epoch: %d\n', epciters);
        fprintf('Learning rate / decay / maxfails: %g %g %g\n', rate, decay, maxfails);
        fprintf('Function Sampler: %s\n', fsampler_str);
        fprintf('Gradient Sampler: %s\n', gsampler_str);
    else
        fprintf('Max iterations: %d\n', maxiters);
        fprintf('Projected gradient tolerance: %.4g\n', pgtol);
    end
    fprintf('\n');
end

%% L-BFGS-B Optimization
if use_lbfgsb

    setupStartB = tic;
    
    fcn = @(x) tt_gcp_fg(update(M0,1:nd,x), X, fh, gh, W, true, true, true);
    
    lbvec = lb * ones(vecsz,1);
    ubvec = inf(vecsz,1);
    
    lbfgsbopts = struct;
    lbfgsbopts.x0 = tovec(M0,false);
    lbfgsbopts.printEvery = printitn;
    lbfgsbopts.maxIts = maxiters;
    lbfgsbopts.maxTotalIts = maxiters*10;
    lbfgsbopts.factr = factr;
    lbfgsbopts.pgtol = pgtol;
    setupTimeB = toc(setupStartB);
    
    if (printitn > 0)
        fprintf('Begin Main loop\n');
    end
    
    mainStart = tic;
    lbfgsbopts.errFcn = @(x) toc(mainStart);
    [x,finalf,lbfgsout] = lbfgsb(fcn, lbvec, ubvec, lbfgsbopts);
    mainTime = toc(mainStart);
    
    M = update(M0,1:nd,x);
       
    info.fcn = fcn;
    info.lbfgsbopts = lbfgsbopts;
    info.lbfgsout = lbfgsout;
    info.finalf = finalf;
    
    if printitn > 0
        fprintf('End Main Loop\n');
        fprintf('\n');
        fprintf('Final objective: %10.4e\n', finalf);
        fprintf('Setup time: %.2f seconds\n', setupTimeA+setupTimeB);
        fprintf('Main loop time: %.2f seconds\n', mainTime);
        fprintf('Outer iterations: %d\n', lbfgsout.iterations);
        fprintf('Total iterations: %d\n', lbfgsout.totalIterations);
        fprintf('L-BFGS-B Exit message: %s\n', lbfgsout.lbfgs_message1);
    end
    
end

%% Stochastic Optimization    
if use_stoc
    
    setupStartB = tic;
    
    % Initialize moments
    if use_adam
        m = cell(nd,1);
        v = cell(nd,1);
        for k = 1:nd
            m{k} = zeros(sz(k),nc);
            v{k} = zeros(sz(k),nc);
        end
    else
        m = [];
        v = [];
    end
    
    % Only used by Adagrad
    gnormsum = 0;
    
    % Extract samples for estimating function value - these never change
    [fsubs,fvals,fwgts] = fsampler();
    
    % Compute initial estimated function value
    fest = tt_gcp_fg_est(M0,fh,gh,fsubs,fvals,fwgts,true,false,false,false);
    
    % Set up loop variables
    M = M0;                         % Copy initial guess
    nfails = 0;                     % Counter # times fails to decrease
    titers = 0;                     % Total iterations (used by ADAM)
    
    Msave = M;                      % Best model so far
    msave = m;                      % Corresponding ADAM parameters
    vsave = v;                      % Corresponding ADAM parameters
    fest_prev = fest;               % Corresponding function value
    
    % Tracing the progress in the function value by epoch
    fest_trace = zeros(maxiters+1,1);
    step_trace = zeros(maxiters+1,1);
    time_trace = zeros(maxiters+1,1);
    fest_trace(1) = fest;

    % Print status
    if (printitn > 0)
        fprintf('Begin Main loop\n');
        fprintf('Initial f-est: %e\n', fest);
    end
    
    setupTimeB = toc(setupStartB);
    mainStart = tic;
    time_trace(1) = toc(setupStartA);
    
    
    % Main loop outer iteration
    for nepoch = 1:maxiters      
        
        % Main loop inner iteration
        step = decay^nfails * rate;     % Ignored by Adagrad
        for iter = 1:epciters
            
            % Tracking total iterations
            titers = titers + 1;
            
            % Select subset for stochastic gradient
            [gsubs,gvals,gwts] = gsampler();
            
            % Compute gradients for each mode
            Gest = tt_gcp_fg_est(M,fh,gh,gsubs,gvals,gwts,false,true,false,false,crng);
            
            % Check for inf
            isinfgrad = cellfun(@(x) any(isinf(x(:))), Gest, 'UniformOutput', true);
            if any(isinfgrad)
                error('Infinite gradient encountered! (epoch = %g, iter = %g)', nepoch, iter);
            end
            
            % Take a step
            if use_adam
                m = cellfun(@(mk,gk) beta1*mk + (1-beta1)*gk,m,Gest,'UniformOutput',false);
                v = cellfun(@(vk,gk) beta2*vk + (1-beta2)*gk.^2,v,Gest,'UniformOutput',false);
                mhat = cellfun(@(mk) mk/(1-beta1^titers),m,'UniformOutput',false);
                vhat = cellfun(@(vk) vk/(1-beta2^titers),v,'UniformOutput',false);
                M.u = cellfun(@(uk,mhk,vhk) max(lb,uk-step*mhk./(sqrt(vhk)+epsilon)),M.u,mhat,vhat,'UniformOutput',false);
            elseif use_adagrad
                gnormsum = gnormsum + sum(cellfun(@(gk) sum(gk(:).^2), Gest, 'UniformOutput',true));
                step = 1/sqrt(gnormsum); 
                M.u = cellfun(@(uk,gk) max(lb,uk-step*gk),M.u,Gest,'UniformOutput',false);
            else
                M.u = cellfun(@(uk,gk) max(lb,uk-step*gk),M.u,Gest,'UniformOutput',false);
            end
        end
        
        % Estimate objective function value
        fest = tt_gcp_fg_est(M,fh,gh,fsubs,fvals,fwgts,true,false,false,false);
        
        % Save trace
        fest_trace(nepoch+1) = fest;
        step_trace(nepoch+1) = step;
        
        % Check convergence condition
        failed_epoch = fest > fest_prev;
        
        if failed_epoch
            nfails = nfails + 1;
        end   
        
        festtol_test = fest < festtol;
        
        % Reporting
        if (printitn > 0) && (mod(nepoch,printitn) == 0 || failed_epoch || festtol_test)
            fprintf('Epoch %2d: f-est = %e, step = %g', nepoch, fest, step);
            if failed_epoch
                fprintf(', nfails = %d (resetting to solution from last epoch)', nfails);
            end           
            fprintf('\n');
            
        end
 
        if failed_epoch
                       
            % Back up to best solution so far!
            M = Msave;
            m = msave;
            v = vsave;
            fest = fest_prev;
            titers = titers - epciters;
            
            % Reset Adagrad
            gnormsum = 0;
            
        else
            
            % Save current solution
            Msave = M;
            msave = m;
            vsave = v;
            fest_prev = fest;
            
        end

        % Save time
        time_trace(nepoch+1) = toc(setupStartA);                
        
        if (nfails > maxfails) || festtol_test
            break;
        end
                
        
    end
    
    mainTime = toc(mainStart);
    
    info.fest_trace = fest_trace(1:nepoch+1);
    info.step_trace = step_trace(1:nepoch+1);
    info.time_trace = time_trace(1:nepoch+1);
    info.nepoch = nepoch;
    
    if printitn > 0
        fprintf('End Main Loop\n');
        fprintf('\n');
        fprintf('Final f-est: %10.4e\n', fest);
        fprintf('Setup time: %.2f seconds\n', setupTimeA+setupTimeB);
        fprintf('Main loop time: %.2f seconds\n', mainTime);
        fprintf('Total iterations: %d\n', nepoch * epciters);
    end
    
    
end

%% Wrap up

% Save timings
info.mainTime = mainTime;
info.setupTimeA = setupTimeA;
info.setupTimeB = setupTimeB;
info.setupTime = setupTimeA+setupTimeB;

% Arrange the final tensor so that the columns are normalized.
M = fixsigns(arrange(M));









