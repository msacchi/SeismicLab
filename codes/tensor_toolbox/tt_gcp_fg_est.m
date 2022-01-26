function [F, G] = tt_gcp_fg_est(M, fh, gh, subs, xvals, weights, computeF, computeG, vectorG, LambdaCheck, crng)
%TT_GCP_FG_EST Estimate the GCP function and gradient with a subsample.
%
%   [F,G] = TT_GCP_FG_EST(M, FH, GH, XSUBS, XVALS, WVALS) estimates the GCP
%   function and gradient specified by FH and GH for M and X. In this case,
%   we have only a portion of X as  specified by XSUBS and XVALS along with
%   the corresponding sampling weights in WVALS that are used in the estimate.
%
%   See also GCP_OPT, TT_GCP_FG, TT_GCP_FG_SETUP.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

% Created by Tamara G. Kolda, Fall 2018. Includes work with
% collaborators David Hong and Jed Duersch. 

%% Hidden options
%
%  Note that there are five hidden options. The first three are similar to
%  the hidden options for gcp_fg and the fourth is whether or not to verify
%  that M has lambda = [1,1,...,1], which is assumed and so should
%  generally be checked unless the user is absolutely sure it's okay. 

%% Parse inputs
if nargin < 11
       
    if ~exist('computeF','var')
        computeF = true;
    end
    
    if ~exist('computeG','var') 
        computeG = (nargout > 1);
    end
    
    if ~exist('vectorG','var')
        vectorG = false;
    end
    
    if ~exist('LambdaCheck','var')
        LambdaCheck = true;
    end
    
    % Specify range for correction/adjustment when nonzeros may be included
    % in the "zero" sample. In this case, crng should be the indices
    % of nonzero samples, which are the ones that are adjusted. 
    if ~exist('idx','var')
        crng = [];
    end
    
end

%% Input checks (keep minimal for timing's sake)

d = ndims(M);
sz = size(M);
F = [];
G = [];

if LambdaCheck && ~all(M.lambda == 1)
    warning('Fixing M to have all ones for lambda');
    M = normalize(M,1);
end

%% Compute model values and exploded Zk matrices
[mvals, Zexp] = gcp_fg_est_helper(M.u, subs);

%% Compute function value
if computeF
    Fvec = fh(xvals,mvals);
    if ~isempty(crng)
        Fvec(crng) = Fvec(crng) - fh(0,mvals(crng));
    end
    F = sum( weights .* Fvec );
end
if ~computeG
    return;
end

%% Compute sample y values
yvals = weights .* gh(xvals, mvals);
if ~isempty(crng)
    yvals(crng) = yvals(crng) - weights(crng) .* gh(0, mvals(crng));
end

%% Compute function and gradient
G = cell(d,1);
nsamples = size(subs,1);
for k=1:d
    % The row of each element is the row index to accumulate in the
    % gradient. The column indices are corresponding samples. They are
    % in order because they match the vector of samples to be
    % multiplied on the right.    
    S = sparse(subs(:,k), (1:nsamples)', yvals, sz(k), nsamples, nsamples);    
    G{k} = S * Zexp{k};
end

% Convert to single vector
if vectorG
    G = cell2mat(cellfun(@(x) x(:), G, 'UniformOutput', false));
end

%% If not computing F, set F (the 1st return arugment) to be the gradient
if ~computeF
    F = G;
end

function [mvals, Zexp] = gcp_fg_est_helper(factors, subs)
% GCP_FG_EST_HELPER Model values at sample locations and exploded Zk's.  

% Created by Tamara G. Kolda, Sept. 2018. Includes prior work by
% collaborators David Hong and Jed Duersch.  

% Check for empty
if isempty(subs)
    mvals = [];
    return;
end

% Process inputs
d = size(subs,2);

% Create exploded U's from the model factor matrices
Uexp = cell(d,1);
for k = 1:d
    Uexp{k} = factors{k}(subs(:,k),:);
end

% After this pass,
% Zexp{k} = Hadarmard product of Uexp{1} through Uexp{k-1}
% for k = 2,...,d.
Zexp = cell(1,d);
Zexp{2} = Uexp{1};
for k = 3:d
    Zexp{k} = Zexp{k-1} .* Uexp{k-1};
end

% After this pass,
% Zexp{k} = Hadamard product of Uexp{1} though Uexp{d}, except Uexp{k}
% for k = 1,...,d.
Zexp{1} = Uexp{d};
for k=d-1:-1:2
    Zexp{k} = Zexp{k} .* Zexp{1};
    Zexp{1} = Zexp{1} .* Uexp{k};
end

% Compute model values at sample locations
mvals = sum(Zexp{d} .* Uexp{d},2);


