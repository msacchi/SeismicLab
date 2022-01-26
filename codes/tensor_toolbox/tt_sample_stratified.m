function [subs, xvals, weights] = tt_sample_stratified(X, xnzidx, nnzs, nzrs, oversample)
%TT_SAMPLE_STRATIFIED Sample nonzero and zero entries from a sparse tensor.
%
%   [SUBS,VALS,WGTS] = TT_SAMPLE_STRATIFIED(X,NZIDX,NNZ,NZR) creates a
%   stratifies sample of nonzero and zero entries of the sparse tensor X.
%   The NZIDX is the sorted 64-bit linear indices of the nonzeros in X.
%   This is required for efficients of the sampling. The values NNZ and NZR
%   specify the desired number of nonzero and zero samples, respectively.
%
%   Example
%      nzidx = tt_sub2ind64(sz,X.subs);
%      nzidx = sort(nzidx);
%      [subs,vals,wgts] = tt_sample_stratified(X,nzidx,1000,1000);
%      [f,G] = tt_gcp_fg_est(M,fh,gh,subs,vals,wgts);
%
%   See also GCP_OPT, TT_GCP_FG_EST, TT_SAMPLE_UNIFORM.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

% Created by Tamara G. Kolda, Fall 2018. Includes work with
% collaborators David Hong and Jed Duersch. 

%% Setup
sz = size(X);
nelx = prod(sz); % Number of elements in X
nnzx = nnz(X); % Number of Nonzeros in X
nzrx = nelx - nnzx; % Number of Zeros in X
with_replacement = true;

if nargin < 5
    oversample = 1.1;
end

%% Sample nonzeros
[nonzero_subs, nonzero_xvals] = tt_sample_nonzeros(X,nnzs,with_replacement); 
nonzero_weights = (nnzx / nnzs) * ones(nnzs,1); 

%% Sample zeros
zero_subs = tt_sample_zeros(X,xnzidx,nzrs,oversample,with_replacement);
zero_xvals = zeros(nzrs,1);
zero_weights = (nzrx / nzrs) * ones(nzrs,1);

%% Assemble "Sample"
subs = [nonzero_subs; zero_subs];
xvals = [nonzero_xvals; zero_xvals];
weights = [nonzero_weights; zero_weights];


