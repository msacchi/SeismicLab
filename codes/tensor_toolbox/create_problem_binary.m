function [X,Mtrue,info] = create_problem_binary(sz,r,varargin)
%CREATE_PROBLEM_BINARY Creates random low-rank 0/1 tensor.
%
%   [X,M,INFO] = CREATE_PROBLEM_BINARY(SZ,R,'param','value') creates an
%   sptensor X of size SZ from the low-rank ktensor M of rank R that
%   corresponds to the *odds* of a 1 in each position. The parameters that
%   control this are as follows:
%
%      'state' - State of random number generator, for reproducing results.
%      'loprob' - Probability of 'noise' one. Default: 0.01.
%      'hiprob' - Probability of 'structural' one. Default: 0.90.
%      'density' - Density of structural entries. Default: 1/r.
%      'verbosity' - Output: 0: None, 1: Minimal (default), 2: Detailed.
%      'spgen' - Avoid explicitly forming low-rank tensor. Default: False.
%
%   REFERENCES: 
%   * T. G. Kolda, D. Hong, J. Duersch. Stochastic Gradients for
%     Large-Scale Tensor Decomposition, 2019.
%
%   See also: GCP_OPT, CREATE_PROBLEM.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

% Created by Tamara G. Kolda, Fall 2018. Includes work with

%% Random set-up
defaultStream = RandStream.getGlobalStream;

%% Set algorithm parameters from input or by using defaults
params = inputParser;
params.addParameter('state', defaultStream.State);
params.addParameter('loprob', 0.01, @(x) isscalar(x) && x > 0 && x < 0.1);
params.addParameter('hiprob', 0.9, @(x) isscalar(x) && x > 0 && x < 1);
params.addParameter('density', []);
params.addParameter('verbosity', 1);
params.addParameter('spgen',false);
params.addParameter('Mtrue',[]);
params.parse(varargin{:});
info.params = params.Results;

%% Initialize random number generator with specified state
defaultStream.State = params.Results.state;

%% Extract parameters
loprob = params.Results.loprob;
hiprob = params.Results.hiprob;
density = params.Results.density;
verbosity = params.Results.verbosity;
spgen = params.Results.spgen;
Mtrue = params.Results.Mtrue;

%% Setup
if verbosity > 0
    fprintf('Creating random problem instance\n');
end

%% Set up for creating factor matrices

% Density specifies the density of high values in the first r-1 columns of
% the factor matrices
if isempty(density)
    density = 1/r; 
end 

% Extract the order of the tensor
d = length(sz);

% Convert the high and low probabilities to the dth root of the
% corresponding odds. 
loval = nthroot(loprob/(1-loprob),d);
hival = nthroot(hiprob/(1-hiprob),d);


%% Populate factor matrices
% The first (r-1) columns of each factor matrix is sparse per the
% specified denstiy. The nonzero values are normal distributed around the
% hival odds ration with a standard deviation of 0.5.
% The last column of each factor matrices is dense but low-valued set to
% the constant loval, corresponding to the general noisyness of binary
% observations. 

if isempty(Mtrue)
    A = cell(d,1);
    for k = 1:d
        if r > 1
            A1v = random('Normal', hival, 0.5, [sz(k),r-1]);
            A1p = rand(sz(k),r-1) < density;
            A1 = max(A1v .* A1p, 0);
        else
            A1 = [];
        end
        A2 = loval * ones(sz(k),1);
        A{k} = [A1,A2];
    end
    Mtrue = ktensor(A); % Correct solution
else
    A = Mtrue.u;
    if verbosity > 0
        fprintf('Using user-specified choice for Mtrue\n');
    end
end
%%

if spgen

    % --- Create all-zero sparse tensor ---
    X = sptensor(sz);
    
    % --- Compute big entries of X, which are expected to be few ---
    
    if verbosity > 1
        fprintf('Generating high probability entries...\n');
    end
       
    % Find possible high values correspond to each component
    subs = [];
    for j = 1:r-1
    
        % Identify nonzeros in each mode
        modeidx = cell(d,1);
        for k = 1:d
            tmp = A{k}(:,j);
            modeidx{k} = find(tmp > 0);
        end
        
        % Count nnzs in each factor
        cnts = cellfun(@length, modeidx);
        
        % Compute total number of entries from this factor
        fcnt = prod(cnts);

        if fcnt > 0
            % Create the subscripts of those entries
            csubs = tt_ind2sub(cnts',(1:fcnt)');
            fsubs = zeros(fcnt,d);
            for k = 1:d
                fsubs(:,k) = modeidx{k}(csubs(:,k));
            end
            subs = [subs; fsubs];
        end
        
    end
    
    subs = unique(subs,'rows');
    nhigh_max = size(subs,1);
    if verbosity > 1
        fprintf('\tmax # high entries = %d\n',nhigh_max);
    end

    if nhigh_max > 0       
        % Compute the probablities at those entries
        Mvals = Mtrue(subs);
        Pvals = Mvals ./ (1 + Mvals);
        Xvals = random('Binomial',1,Pvals);
        tf = (Xvals == 1);
        % Remove the subscripts that don't correspond to ones
        hisubs = subs(tf,:);
        X(hisubs) = 1;
        nhigh = sum(Xvals);
    else 
        hisubs = [];
        nhigh = 0;
    end
    
    if verbosity > 1
        fprintf('\t# high entries = %d\n',nhigh);
    end
     
    % --- Compute the 'noise' from the rest of the entries ---
    
    if verbosity > 1
        fprintf('Generating low probability (aka noise) entries...\n');
    end
    
   % Number of remaining entries
    nloprob = prod(sz) - nhigh_max; 
    % Randomly compute how many will be 1's, using binomial,
    % which we estimate using Poisson since nloprob is large and loprob is
    % small.
    nlow = random('Poisson', nloprob * loprob, 1);
    if verbosity > 1
        fprintf('\t# low entries = %d\n',nlow);
    end
    if nlow > 0
        % Choose that many indicies
        losubs = tt_sample_zeros(X, tt_sub2ind(sz,hisubs), nlow, 1.1, false);
        X(losubs) = 1;
    end
    if verbosity > 1
        fprintf('\tFinished\n');
    end

    info.nlow = nlow;
    info.nhigh_max = nhigh_max;
    info.nhigh = nhigh;
    

else
    Mtruef = full(Mtrue);
    P = Mtruef ./ (1 + Mtruef);
    X = sptensor(random('Binomial',1,double(P)));
end
