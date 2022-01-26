function [P, P0, output] = cp_opt(Z,R,varargin)
%CP_OPT Fits a CP model to a tensor via optimization.
%
%   K = CP_OPT(X,R) fits an R-component CANDECOMP/PARAFAC (CP) model
%   to the tensor X. The result K is a ktensor. The function being
%   optimized is F(K) = 1/2 || X - K ||^2.
%
%   K = CP_OPT(X,R,'param',value,...) specifies additional
%   parameters for the method. Specifically...
%
%   'init' - Initialization for factor matrices (default: 'randn'). This
%   can be a cell array with the initial matrices, a ktensor, or one of the
%   following strings:
%      'randn'  Randomly generated via randn function
%      'rand'   Randomly generated via rand function
%      'zeros'  All zeros
%      'nvecs'  Selected as leading left singular vectors of X(n)
%
%   'opt_options' - Optimization method options, passed as a structure.
%   Type 'help lbfgsb' to see the options. (Note that the 'opts.x0' option
%   is overwritten using the choice for 'init', above.) 
%
%   'lower'/'upper' - Lower/upper bounds, passed in as a scalar (if they
%   are all the same), vector, cell array, or ktensor (lambda values
%   ignored).  
%
%   [K, U0] = CP_OPT(...) also returns the initial guess.
%
%   [K, U0, OUT] = CP_OPT(...) also returns a structure with the
%   optimization exit flag, the final relative fit, and the full
%   output from the optimization method. The fit is defined as 
%
%      FIT = 100 * (1 - ( F(K) / F(0) )).
%
%   REFERENCE: E. Acar, D. M. Dunlavy and T. G. Kolda, A Scalable
%   Optimization Approach for Fitting Canonical Tensor Decompositions,
%   J. Chemometrics, 25(2):67-86, 2011, http://doi.org/10.1002/cem.1335.
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','cp_opt_doc.html')))">Documentation page for CP-OPT</a>
%
%   See also TENSOR, SPTENSOR, KTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%% Error checking
% if ~isa(Z,'tensor') && ~isa(Z,'sptensor')
%     error('Z must be a tensor or a sptensor');
% end
% 
if (nargin < 2)
    error('Error: invalid input arguments');
end

%% Set parameters
params = inputParser;
params.addParameter('opt', 'lbfgsb', @(x) ismember(x,{'ncg','tn','lbfgs','lbfgsb'}));
params.addParameter('init', 'randn', @(x) (iscell(x) || isa(x, 'ktensor') || ismember(x,{'random','rand','randn','nvecs','zeros'})));
params.addParameter('lower',-Inf);
params.addParameter('upper',Inf);
params.addParameter('opt_options', '', @isstruct);
params.parse(varargin{:});

init = params.Results.init;
opt = params.Results.opt;
options = params.Results.opt_options;
lower = params.Results.lower;
upper = params.Results.upper;

use_lbfgsb = strcmp(opt,'lbfgsb');

%% Initialization
sz = size(Z);
N = length(sz);

if iscell(init)
    P0 = init;
elseif isa(init,'ktensor')
    P0 = tocell(init);
else
    P0 = cell(N,1);
    if strcmpi(init,'nvecs')
        for n=1:N
            P0{n} = nvecs(Z,n,R);
        end
    else
        for n=1:N
            P0{n} = matrandnorm(feval(init,sz(n),R));
        end
    end
end

%% Set up lower and upper (L-BFGS-B only)

if ~use_lbfgsb && ( any(isfinite(lower)) || any(isfinite(upper)) )
    error('Cannot use lower and upper bounds without L-BFGS-B');
end

if use_lbfgsb
    lower = convert_bound(lower,sz,R);
    upper = convert_bound(upper,sz,R);
end


%% Set up optimization algorithm

if use_lbfgsb % L-BFGS-B
    if ~exist('lbfgsb','file')
        error(['CP_OPT requires L-BFGS-B function. This can be downloaded'...
            'at https://github.com/stephenbeckr/L-BFGS-B-C']);
    end
else % POBLANO
    switch (params.Results.opt)
        case 'ncg'
            fhandle = @ncg;
        case 'tn'
            fhandle = @tn;
        case 'lbfgs'
            fhandle = @lbfgs;
    end
    
    if ~exist('poblano_params','file')
        error(['CP_OPT requires Poblano Toolbox for Matlab. This can be ' ...
            'downloaded at http://software.sandia.gov/trac/poblano.']);
    end     
end

%% Set up optimization algorithm options
if isempty(options) 
    if use_lbfgsb
        options.maxIts = 10000;
        options.maxTotalIts = 50000;
        options.printEvery = 10;
    else
        options = feval(fhandle, 'defaults');
    end
end
    



%% Fit CP using CPOPT
normsqr = norm(Z)^2;
if use_lbfgsb
    opts = options;
    opts.x0 = tt_fac_to_vec(P0);    
    [xx,ff,out] = lbfgsb(@(x)tt_cp_fun(x,Z,normsqr), lower, upper, opts);
    P = ktensor(tt_cp_vec_to_fac(xx, Z));
    output.ExitMsg = out.lbfgs_message1;
    output.Fit = 100 * (1 - ff /(0.5 * normsqr));
    output.OptOut = out;
else % POBLANO
    out = feval(fhandle, @(x)tt_cp_fun(x,Z,normsqr), tt_fac_to_vec(P0), options);
    P = ktensor(tt_cp_vec_to_fac(out.X, Z));
    output.ExitFlag  = out.ExitFlag;
    output.Fit = 100 * (1 - out.F /(0.5 * normsqr));
    output.OptOut = out;
end


%% Clean up final result
% Arrange the final tensor so that the columns are normalized.
P = arrange(P);
% Fix the signs
P = fixsigns(P);


function newbound = convert_bound(bound,sz,R)

len = sum(sz)*R;

if isscalar(bound)
    newbound = bound * ones(len,1);
elseif isa(bound,'ktensor')
    newbound = tovec(bound, false);
elseif iscell(bound)
    newbound = tt_fac_to_vec(bound);
end

if ~isequal(size(newbound), [len 1])
    error('Bound is the wrong size');
end
