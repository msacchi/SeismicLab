function [F,G] = tt_gcp_fg(M, X, f, g, W, computeF, computeG, vectorG)
%TT_GCP_FG Loss function and gradient for generalized CP.
%
%   F = TT_GCP_FG(M,X,f) expects that M is a ktensor and X is a
%   dense tensor. The f is a function handle of the form f(x,m) that
%   measures the elementwise loss for a data entry x and corresponding
%   model entry m. The function should be able to accept vector inputs to
%   do evaluations in bulk.  
%
%   [F,G] = TT_GCP_FG(M,X,f,g) also computes the gradient. Here g is the
%   gradient of f and has the same form. The G is returned as a cell array
%   where G{k} is a matrix that is the same size as M.u{k} (the k-th factor
%   matrix). 
%
%   [F,G] = TT_GCP_FG(M,X,f,g,W) specifies a weight tensor where W is a tensor
%   that is the same size as X. It should have 1's for known values and 0's
%   for missing values. The function/gradient is only computed w.r.t. the
%   known values. Setting W to [] indicated no missing data.
%
%   See also GCP_OPT, TT_GCP_FG_SETUP.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

% Created by Tamara G. Kolda, Fall 2018. Includes work with
% collaborators David Hong and Jed Duersch. 

%% Hidden options
%
%   G = GCP_FG(M,X,f,g,W,false,true) computes only the gradient.
%
%   G = GCP_FG(M,X,f,g,W,false,true,true) computes only the gradient and
%   converts it to a vector (equivalent to the tovec operation on a
%   ktensor).
%
%   [F,G] = GCP_FG(M,X,f,g,W,true,true,true) computes boths the function
%   and the gradient and converts the gradient to vector form.
%

%% Parse inputs
if nargin < 8  
    
    if ~exist('W','var')
        W = [];
    end
    
    if ~exist('computeF','var')
        computeF = true;
    end
    
    if ~exist('computeG','var') 
        computeG = (nargout > 1);
    end
    
    if ~exist('vectorG','var')
        vectorG = false;
    end
    
end

%% Setup
Mfull = full(M);
Mv = Mfull(:);
Xfull = full(X);
Xv = Xfull(:);
F = [];
G = [];

%% Calculate function value
if computeF
    
    Fvec = f(Xv, Mv); % F is a vector
    
    if ~isempty(W)
        Fvec = W(:).*Fvec;  % be sure to zero out any unknown entries
    end
    
    F = sum(Fvec);
    
end

%% QUIT IF ONLY NEED FUNCTION EVAL
if ~computeG
    return;
end

%% Gradient calculation
Y = g(Xv,Mv); % Result is a vector
Y = tensor(Y,size(X));
if ~isempty(W)
    Y = W.*Y;
end

%% Gradient wrt U's using MTTKRP sequence.
G = mttkrps(Y,M.u);

%% Assemble gradient 
if vectorG
    G = cell2mat(cellfun(@(x) x(:), G, 'UniformOutput', false));
end

%% If not computing F, set F (the 1st return arugment) to be the gradient
if ~computeF
    F = G;
end
