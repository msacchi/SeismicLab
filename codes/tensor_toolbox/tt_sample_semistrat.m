function [subs, xvals, weights] = tt_sample_semistrat(X, nnzs, nzrs)
%TT_SAMPLE_SEMISTRAT Sample nonzero and zero entries from a sparse tensor.
%
%   [SUBS,VALS,WGTS] = TT_SAMPLE_OVERLAPPED(X,NNZ,NZR) creates a
%   stratified sample of nonzero and zero entries of the sparse tensor X.
%
%   Example
%   [subs,vals,wgts] = tt_sample_overlapped(X,1000,1000);
%   [f,G] = tt_gcp_fg_est(M,fh,gh,subs,vals,wgts,true,true,true,false,1000);
%
%   See also GCP_OPT, TT_GCP_FG_EST, TT_SAMPLE_UNIFORM.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

% Created by Tamara G. Kolda, Fall 2018. Includes work with
% collaborators David Hong and Jed Duersch. 

%% Setup
d = ndims(X);
sz = size(X);
nelx = prod(sz); % Number of elements in X
nnzx = nnz(X); % Number of Nonzeros in X
nzrx = nelx - nnzx; % Number of Zeros in X
with_replacement = true;

%% Sample nonzeros
[nonzero_subs, nonzero_xvals] = tt_sample_nonzeros(X,nnzs,with_replacement); 
nonzero_weights = (nnzx / nnzs) * ones(nnzs,1); 


%% Sample 'zeros'
zero_subs = bsxfun(@(a,b)ceil(a.*b), rand(nzrs,d), sz);
zero_xvals = zeros(nzrs,1);
zero_weights = nelx / nzrs * ones(nzrs,1);

%% Assemble "Sample"
subs = [nonzero_subs; zero_subs];
xvals = [nonzero_xvals; zero_xvals];
weights = [nonzero_weights; zero_weights];


