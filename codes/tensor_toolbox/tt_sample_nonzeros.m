function [subs,vals] = tt_sample_nonzeros(X, nns, with_replacement)
%TT_SAMPLE_NONZEROS Sample nonzeros from a sparse tensor.
%
%   [VALS,SUBS] = TT_SAMPLE_NONZEROS(X,N) finds N random nonzero entries
%   (uniformly with replacement) in the sparse tensor X. It returns
%   VALS, the values, and SUBS, the corresponding subscripts.  Throws an
%   error if N > nnz(X).
%
%   See also SAMPLE_STRATIFIED.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

% Created by Tamara G. Kolda, Fall 2018. Includes work with
% collaborators David Hong and Jed Duersch. 

% Created by Tamara G. Kolda, Sept. 2018.  

%% Input checks
if nargin < 3
    with_replacement = true;
end

%% Setup
nnx = nnz(X);

%% Select nonzeros
if nns == nnx
    nidx = 1:nnx;
elseif with_replacement
    nidx = randi(nnx, nns,1);
else    
    if nns > nnx
        error('Tensor does not have enough nonzeros to sample');
    end
    nidx = randperm(nnx,nns);
end

%% Extract subscripts and values
subs = X.subs(nidx,:); 
vals = X.vals(nidx); 
