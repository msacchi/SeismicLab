function [fh,gh,lowerbnd] = tt_gcp_fg_setup(type,X)
%TT_GCP_FG_SETUP Sets the GCP functions according to specified name.
%
%   [F,G,L] = TT_GCP_FG_SETUP(TYPE,X) returns the function and gradient
%   function as well as the lower bound for different types of objective
%   functions. It also checks that X satisfies the standard constraints for
%   that choice. The valid types are:
%
%      - 'normal' or 'Gaussian'
%      - 'binary' or 'Bernoulli-odds'
%      - 'Bernoulli-logit'
%      - 'count' or 'Poisson'
%      - 'Poisson-log'
%      - 'Rayleigh'
%      - 'Gamma'
%
%   [F,G,L] = TT_GCP_FG_SETUP(TYPE,X) works for types that require a
%   parameter, which is specified inside the type string. 
%
%      - 'negative-binomial (number of failures)'
%      - 'Huber (delta threshold)'
%      - 'beta-divergence (beta)'
%
%   Details of the functions can be found in D. Hong, T. G. Kolda, J. A. Duersch.
%   Generalized Canonical Polyadic Tensor Decomposition. arXiv:1808.07452,
%   2018.  
%
%   See also GCP_OPT, TT_GCP_FG, TT_GCP_FG_EST.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


tmp = split(type);
type = tmp{1};
if length(tmp)>1
    param = str2double(tmp{2}(2:end-1));
end

switch lower(type)
    case {'normal','gaussian'}
        fh = @(x,m) (m-x).^2;
        gh = @(x,m) 2.*(m-x);
        lowerbnd = -Inf;
    case {'binary', 'bernoulli-odds'}
        if exist('X','var') && ~valid_binary(X)
            warning('Using ''%s'' type but tensor X is not binary', type);
        end
        fh = @(x,m) log(m+1) - x.*log(m + 1e-10);
        gh = @(x,m) 1./(m+1) - x./(m + 1e-10);
        lowerbnd = 0;
    case {'bernoulli-logit'}
        if exist('X','var') && ~valid_binary(X)
            warning('Using ''%s'' type but tensor X is not binary', type);
        end
        fh = @(x,m) log(exp(m) + 1) - x .* m;
        gh = @(x,m) exp(m)./(exp(m) + 1) - x;
        lowerbnd = -Inf;
    case {'count','poisson'}
        if exist('X','var') && ~valid_natural(X)
            warning('Using ''%s'' type but tensor X is not counts', type);
        end
        fh = @(x,m) m - x.*log(m + 1e-10);
        gh = @(x,m) 1 - x./(m + 1e-10);
        lowerbnd = 0;
    case {'poisson-log'}
         if exist('X','var') && ~valid_natural(X)
            warning('Using ''%s'' type but tensor X is not counts', type);
        end
        fh = @(x,m) exp(m) - x.*m;
        gh = @(x,m) exp(m) - x;
        lowerbnd = -Inf;        
    case 'rayleigh'
        if exist('X','var') && ~valid_nonneg(X)
            warning('Using ''%s'' type but tensor X is not nonnegative', type);
        end
        fh = @(x,m) 2*log(m+1e-10) + (pi/4)*(x./(m+1e-10)).^2;
        gh = @(x,m) 2./(m+1e-10) - (pi/2)*x.^2./(m+1e-10).^3;
        lowerbnd = 0;
    case 'gamma'
        if exist('X','var') && ~valid_nonneg(X)
            warning('Using ''%s'' type but tensor X is not nonnegative', type);
        end
        fh = @(x,m) x./(m+1e-10) + log(m+1e-10);
        gh = @(x,m) -x./((m+1e-10).^2) + 1./(m+1e-10);
        lowerbnd = 0;
    case 'huber'
        if ~exist('param','var')
            error('Need to specify threshold')
        end
        d = param;
        eval(sprintf('fh = @(x,m) (x-m).^2  .* (abs(x-m) < %g) + (%g .* abs(x-m)- %g) .* (abs(x-m) >= %g);',d,2*d,d^2,d));
        eval(sprintf('gh = @(x,m) -2.*(x-m) .* (abs(x-m) < %g) - (%g.*sign(x-m))  .* (abs(x-m) >= %g);',d,2*d,d));
        lowerbnd = -Inf;
    case 'negative-binomial'
        if exist('X','var') && ~valid_nonneg(X)
            warning('Using ''%s'' type but tensor X is not nonnegative', type);
        end
        if ~exist('param','var')
            error('Need to specify number of trials')
        end
        r = param;
        eval(sprintf('fh = @(x,m) (%d+x) .* log(1+m) - x * log(m+1e-10);',r));
        eval(sprintf('gh = @(x,m) (%d)./(1+m) - x./(m+1e-10);',r+1));
        lowerbnd = 0;
     case 'beta'
        if exist('X','var') && ~valid_nonneg(X)
            warning('Using ''%s'' type but tensor X is not nonnegative', type);
        end
        if ~exist('param','var')
            error('Need to specify beta')
        end
        b = param;
%        eval(sprintf('fh = @(x,m) (1/%g) .* (m+1e-10).^%g - (1/(%g-1)) .* x .* (m+1e-10).^(%g-1);',b,b,b,b));
%        eval(sprintf('gh = @(x,m) (m+1e-10).^(%g-1) - x.*(m+1e-10).^(%g-2);',b,b));
        eval(sprintf('fh = @(x,m) (%g) .* (m+1e-10).^(%g) - (%g) .* x .* (m+1e-10).^(%g);',1/b,b,1/(b-1),b-1));
        eval(sprintf('gh = @(x,m) (m+1e-10).^(%g) - x.*(m+1e-10).^(%g);',b-1,b-2));
        lowerbnd = 0;
    otherwise
        error('Unknown type: %s', type);
end

function tf = valid_nonneg(X)

if isa(X,'sptensor')
    tf = all(X.vals > 0);
else
    tf = all(X(:) > 0);
end

function tf = valid_binary(X)

if isa(X,'sptensor')
    tf = all(X.vals == 1);
else
    tf = isequal(unique(X(:)),[0;1]);
end

function tf = valid_natural(X)

if isa(X, 'sptensor')
    vals = X.vals;
else
    vals = X(:);
end

tf = all(vals >= 0) && all(vals == round(vals));
