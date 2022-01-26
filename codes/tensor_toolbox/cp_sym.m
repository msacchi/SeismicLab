function [Model,Info] = cp_sym(A,P,varargin)
%CP_SYM Fit a symmetric P model to the symmetric input tensor.
%   
%   MODEL = CP_SYM(X,R) fits an R-component symmetric CP model to the
%   tensor X. The result MODEL is a symmetric Kruskal tensor (symktensor).
%   The function being optimized is defined by the parameters listed
%   below.
%
%   MODEL = CP_SYM(X,R,'param',value) accepts additional arguments as
%   described below.
%
%   Several optimization methods can be used, depending on what toolboxes
%   and codes are available. The MATLAB optimization toolbox provides
%   FMINUNC and FMINCON. The POBLANO package provides limited-memory BFGS
%   (LBFGS), Nonlinear CG (NCG), and Truncated Newton (TN).
%
%   o 'alg'         - Optimiation algorithm. Choices: {'fminunc','fmincon',
%                     'lbfgs','ncg','tn'}. Default: 'lbfgs'.
%   o 'alg_options' - Options that are passed to the optimization
%                     algorithm. 
% 
%   Parameters that define the objective function (more details in
%   symktensor/setup_fg)... 
%
%   o 'unique'    - Give each unique index equal weight. Default: True.
%   o 'fast'      - Use fast version if unique is false. Default: True.
%   o 'l2weight'  - Penalty on column norms of X = 1. Default: 0.  
%   o 'l1weight'  - Weight to encourage sparsity in LAMBDA. Default: 0.
%   o 'l1param'   - Parameter to encourage sparsity in LAMBDA. Default: 10.
%   o 'nonneg'    - Require: X >= 0, LAMBDA >= 0. Default: False.
%   o 'nolambda'  - Remove LAMBDA from the optimization. Default: False.
%
%   Additionally, it is possible to specify an initial guess.
% 
%   o 'init'      - Specify initial guess. Default: [] (none).
%
%   [MODEL,INFO] = CP_SYM(X,R,...) returns additional information about the
%   optimization: 
%
%   o model0 - Initial guess for model.
%   o data - Produced by fg_setup using X and the initial guess.
%   o setuptime - Time for setup.
%   o optout - Information returned by the optimization method.
%   o optalg - Optimization algorithm. (See 'alg' above.)
%   o optopt - Optimization parameters. (See 'alg_options' above.)
%   o runtime - Time for running optimization method.
%
%   Reference: T. G. Kolda, Numerical Optimization for Symmetric Tensor
%   Decomposition, Mathematical Programming B, Vol. 151, No. 1, pp.
%   225-248, 2015, https://doi.org/10.1007/s10107-015-0895-0 
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','cp_sym_doc.html')))">Documentation page for CP-SYM</a>
%
%   See also SYMKTENSOR, TENSOR/ISSYMMETRIC, SYMKTENSOR/FG,
%   SYMKTENSOR/FG_SETUP.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

%% Check inputs
if nargin < 2
    error('Requires at least two arguments')
end
if ~isa(A,'tensor') && ~isa(A,'symtensor')
    error('First input must be a tensor or symtensor');
end
if ~issymmetric(A)
    error('Input tensor must be symmetric');
end
if ~isscalar(P) || (P < 1) || (round(P) ~= P)
    error('Second input must be a scalar positive integer value');
end

%% Parse optional parameters
params = inputParser;
params.addParameter('alg','lbfgs',@(x) ismember(x,{'fminunc','fmincon','lbfgs','ncg','tn','snopt'}));
params.addParameter('alg_options',[]);
params.addParameter('init',[]);
params.KeepUnmatched = true; % <-- Unmatched params passed to fg_setup
params.parse(varargin{:});

alg = params.Results.alg;
options = params.Results.alg_options;
Model0 = params.Results.init;

%% Setup
if isempty(Model0)
    Model0 = symktensor(P,A); %random initialization. With rand, i.e. U(0,1)
elseif isnumeric(Model0)
    Model0 = symktensor(ones(P,1),Model0,P);
end

starttime = tic;
data = fg_setup(Model0,A,params.Unmatched);
x0 = tovec(Model0,data.fgopts.nolambda);
fhandle = @(x) fg_wrapper(x,data);
setuptime = toc(starttime);

Info.model0 = Model0;
Info.data = data;
Info.setuptime = setuptime;

%% Optimization
if isequal(alg(1:4),'fmin') % MATLAB Optimization Toolbox

    if isempty(options)
        options = optimset('GradObj','on',...
            'MaxFunEvals',50000,...
            'MaxIter',10000,...
            'TolFun',1e-8,...
            'Display','iter');
    end
        
    if  data.fgopts.nonneg
        lb = zeros(size(x0));
        ub = Inf * ones(size(x0));
        starttime = tic;
        [x,fval,flag,out] = fmincon(fhandle,x0,[],[],[],[],lb,ub,[],options);
    else
        starttime = tic;
        [x,fval,flag,out] = fminunc(fhandle,x0,options);
    end
    
    runtime = toc(starttime);    
    Model = symktensor(x,A,data.fgopts.nolambda);
    out.fval = fval;
    out.flag = flag;
    Info.optout = out;
    Info.optalg = alg;
    Info.optopt = options;    
    Info.runtime = runtime;

elseif ismember(alg,{'lbfgs','ncg','tn'}) % POBLANO

    if (data.fgopts.nonneg)
        error('Nonnegative optimization requires a different method');
    end
    optfh = eval(sprintf('@%s',alg));
    if isempty(options)
        options = optfh('defaults');
        options.DisplayIters = 10;
        options.Display = 'iter';
        options.StopTol = 1e-7;
        options.RelFuncTol = 1e-8;
        options.MaxIters = 10000;
        options.MaxFuncEvals = 50000;
    end  
    starttime = tic;
    out = optfh(fhandle, x0, options);
    runtime = toc(starttime);
    Model = symktensor(out.X,A,data.fgopts.nolambda);
    Info.optout = out;
    Info.optalg = alg;
    Info.optopt = options;    
    Info.runtime = runtime;

elseif ismember(alg,{'snopt'}) % SNOPT

    % NOTE that the MATLAB interface for SNOPT requires that we create (and
    % delete) some temporary files. Specifically...
    %
    %  o snoptspecs.txt 
    %  o snoptoutput.txt 
    %  o snoptwrapper.m
    %
    % Also, the way the options works are different for SNOPT than the
    % other methods. The user can define an options struct with the
    % parameter names for SNOPT, replacing spaces with underscores.
    
    if isempty(options)
        options.Major_Iteration_limit = 10000;
        options.New_superbasics_limit = 999;
        options.Superbasics_limit = 999;
        options.Major_optimality_tolerance = 1e-8;
    end
    
    specfile = which('snoptspecs.txt');
    wd_specfile = isempty(specfile);
    if wd_specfile
        fid = fopen('snoptspecs.txt','w');
        fprintf(fid,'Begin snmain2\n');
        fprintf(fid,'   Derivative option                1\n');
        fprintf(fid,'   Major iterations               200\n');
        fprintf(fid,'   Major Print level           000001\n');
        fprintf(fid,'*                             (JFLXBT)\n');
        fprintf(fid,'   Minor print level                1\n');
        fprintf(fid,'   Solution                        No\n');
        fprintf(fid,'End snmain2\n');
        fclose(fid);
        specfile = which('snoptspecs.txt');
    else
        warning('Using existing ''snoptspecs.txt'' file.');
    end
    if isempty(specfile)
        error('SNOPT specification file not found');
    end
    if exist('snoptoutput.txt','file')
        warning('Overwriting ''snoptoutput.txt'' file.');
    end
    snprint('snoptoutput.txt');
    snspec ( specfile );    
    names = fieldnames(options);
    for i = 1:length(names)
        n = names{i};
        estr = sprintf('snseti (''%s'', %g);', strrep(n,'_',' '), options.(n));
        eval(estr);
    end
    snset  ('Minimize');
    wd_wrapper = ~exist('snoptwrapper.m','file');
    if wd_wrapper
        fid = fopen('snoptwrapper.m','w');
        fprintf(fid,'function [F,G] = snoptwrapper(x)\n');
        fprintf(fid,'global snoptwrapperfun\n');
        fprintf(fid,'[F,G] = snoptwrapperfun(x);\n');
        fclose(fid);
    else
        warning('Using existing ''snoptwrapper.m'' file.');
    end
    
    global snoptwrapperfun
    snoptwrapperfun = @(x) fg_wrapper(x,data);
    
    lb = -Inf * ones(size(x0));
    ub = Inf * ones(size(x0));
    if  data.fgopts.nonneg
        lb = zeros(size(x0));
    end
    
    starttime = tic;
    [x,fval,flag] = snopt(x0,lb,ub,0,Inf,'snoptwrapper');
    snprint off; % Closes the file and empties the print buffer
    runtime = toc(starttime);
    
    %clear global fhandle
    Model = symktensor(x,A,data.fgopts.nolambda);
    out.fval = fval;
    out.flag = flag;
    Info.optout = out;
    Info.optoptions = options;
    Info.optalg = alg;
    Info.runtime = runtime;
    
    delete('snoptoutput.txt');
    if wd_specfile
        delete('snoptspecs.txt');
    end
    if wd_wrapper
        delete('snoptwrapper.m');
    end
else
    
    error('Invalid optimization method');
    
end


%% Extract solution and Info about the run

function [f,g] = fg_wrapper(x,data)
model = symktensor(x,data.M,data.P,data.fgopts.nolambda);
[f,g] = fg(model,data);
