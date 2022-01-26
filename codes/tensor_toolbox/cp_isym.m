function [M,info,M0,varargin] = cp_isym(X,r,varargin)
%CP_ISYM Implicit symmetric CP for tensors formed from outer products.
%
%   M = CP_ISYM(X,R) requires X to be a symktensor and R is the number
%   of desired components in the resulting symmetric CP stored in the
%   output symktensor M. Typically, R << ncomponents(X). 
%
%   [M,INFO] = CP_ISYM_OPT(X,R) returns additional information in INFO.
%
%   [...] = CP_ISYM_OPT(X,R,'param','value') takes additional arguments:
%
%      'method' - Optimzation algorithm. Default:'lbfgsb'. 
%         o 'lbfgsb' (Quasi-Newton method with bound constraints),
%         o 'lbfgs' (Quasi-Newton method from Poblano Toolbox), 
%         o 'fminunc'(Quasi-Netwon from Optimization Toolbox),
%         o 'adam' (stochastic gradient descent with momentum).              
%         For optimzation algorithm choices and parameters, see
%         <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','opt_options_doc.html')))">Documentation for Tensor Toolbox Optimization Methods</a>
%
%      'init' - Initialization for factor matrix, Default: 'rrf'. 
%         The user can pass in an explicit initial guess as a symktenor or
%         a factor matrix, pass in a function handle to generate the matrix
%         given the size as two arugments, or specify 'rrf' for the
%         randomized range finder.  
%
%      'state' - Random state, to re-create the same outcome.
%
%      'Xnormsqr' - Objective function f(M) = ||X||^2 + ||M||^2 - 2 <X,M>
%         has a constant term (||X||^2) that can be ignored.
%         Choices are 'exact'=compute or specified value. Default: 0. 
%
%   Reference: S. Sherman, T. G. Kolda, Estimating Higher-Order Moments
%   Using Symmetric Tensor Decomposition, SIAM J. Matrix Analysis and
%   Applications, 41:1369-1387, 2020, https://doi.org/10.1137/19m1299633 
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','cp_isym_doc.html')))">Additional documentation for CP-ISYM</a>
%
%   See also SYMKTENSOR/FG_IMPLICIT.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

%%
tic; % Start timer for setup costs - finished inside optimization call!

%% Error Checking
if ~isa(X,'symktensor')
    error('Ximp must be a symktensor');
end

%% Random set-up
defaultStream = RandStream.getGlobalStream;

%% Algorithm Parameters
params = inputParser;
params.KeepUnmatched = true;
params.PartialMatching = false;
params.addParameter('state', defaultStream.State);
params.addParameter('init','rrf');
params.addParameter('method','lbfgsb');
params.addParameter('Xnormsqr', 0);
params.addParameter('fsamp', 1000);
params.addParameter('gsamp', 100);
params.addParameter('printitn',1);

params.parse(varargin{:});

%% Initialize random number generator with specified state
defaultStream.State = params.Results.state;

%% Setup
init = params.Results.init;
method = params.Results.method;
Xnormsqr = params.Results.Xnormsqr;
fsamp = params.Results.fsamp;
gsamp = params.Results.gsamp;
printitn = params.Results.printitn;

optopts = params.Unmatched;
optopts.printitn = printitn;

% Save
info.params = params.Results;

%% Extract sizes, etc.
nd = ndims(X);
n = size(X,1);

%% Initial guess
if isa(init,'symktensor')
    M0 = init;
    info.params.init = '[deleted]';
else
    if isa(init,'function_handle')
        A0 = init(n,r);
    elseif strcmpi(init,'rrf')
        A0 = matrandnorm( X.u * randn(size(X.u, 2), r) );
    elseif isequal(size(init),[n,r])
        A0 = init;
        info.params.init = '[deleted]';
    else
        error('Invalid ''init'' option');
    end
    lambda0 = ones(r,1);
    M0 = symktensor(lambda0,A0,nd);
end
if ncomponents(M0) ~= r
    error('Initial guess has %d components but expected %d components', ncomponents(M0), r);
end

%% Compute ||X||^2
if ~isscalar(Xnormsqr)
    if isempty(Xnormsqr)
        Xnormsqr = 0;
    elseif strcmpi(Xnormsqr,'exact')
        Xnormsqr = norm(X)^2;
    end
end

%% Shared options
optopts.xdesc = sprintf('Order, Size: %d, %d', ndims(X), size(X,1));

%% Finish setup
setuptime = toc; 


%% Optimization
tic
if (printitn > 0)
    fprintf('\nCP-ISYM Implict Symmetric CP Optimzation');
end
switch(method)
    case {'lbfgsb','lbfgs','fminunc'}
        fgh = @(x) fg_implicit(update(M0,x),X,Xnormsqr);
        optname = sprintf('tt_opt_%s',method);
        [x,f,optinfo] = feval(optname, tovec(M0), fgh, optopts);

    case {'adam'}
        if strcmpi(fsamp,'exact')
            optopts.fexact = true;
            XF = X;
            optopts.fdesc = 'Function: exact';
        else
            optopts.fexact = false;
            XF = randextract(X,fsamp);
            if isscalar(fsamp)
                optopts.fdesc = sprintf('Function: %d samples out of %d observations', fsamp, ncomponents(X));
            else
                optopts.fdesc = sprintf('Function: %d user-specified samples out of %d observations', length(fsamp), ncomponents(X));
            end
        end
        fh = @(mvec) f_implicit(update(M0,mvec),XF,Xnormsqr);
        gh = @(mvec) g_implicit(update(M0,mvec), randextract(X,gsamp));
        optopts.gdesc = sprintf('Gradient: using %d samples  out of %d observations', gsamp, ncomponents(X));       
        [x,f,optinfo] = tt_opt_adam(tovec(M0), fh, gh, optopts);
        
    otherwise
        error('Invalid method')
end
opttime = toc;
%% Clean up
M = normalize(update(M0,x));

%% Save results

info.f = f;
info.optout = optinfo;
info.opttime = opttime;
info.setuptime = setuptime;
info.Xnormsqr = Xnormsqr;
info.opttime = opttime;
info.setuptime = setuptime;





