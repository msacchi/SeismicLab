function subs = tt_sample_zeros(X,Xnzidx,nsamps,oversample,with_replacement)
%TT_SAMPLE_ZEROS Sample zero entries from a sparse tensor.
%
%   SUBS = TT_SAMPLE_ZEROS(X,NZIDX,N) finds N random zero entries (uniformly
%   with replacement) in the sparse tensor X, where the zeros are not
%   stored explicitly. The NZIDX is the sorted linear indices of the
%   nonzeros in X. The return value SUBS is a list of subscripts of zero
%   entries. The procedure automatically determines how much it needs to 
%   oversample to find N nonzero entries and prints a warning if it fails
%   to do so.
%
%   SUBS = TT_SAMPLE_ZEROS(X,N,OR) is the same as above except that OR > 1
%   specifies the oversample rate. The default is 1.1, but this can be
%   increased if the method has trouble getting enough samples which can
%   happen when N is relatively large compared ot the number of zeros in X.
%
%   See also SAMPLE_STRATIFIED.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

% Created by Tamara G. Kolda, Fall 2018. Includes work with
% collaborators David Hong and Jed Duersch. 
 

%% Setup
d = ndims(X);
sz = size(X);
nelx = prod(sz); % Number of entries in X
nnzx = length(Xnzidx); % Number of Nonzeros in X
nzrx = nelx - nnzx; % Number of Zeros in X

if nargin < 3
    oversample = 1.1;
elseif oversample < 1.1
    error('Oversampling rate must be >= 1.1');
end

if nargin < 4
   with_replacement = true;    
end


%% Determine the number of samples to generate
% We need to oversample to account for potential duplicates and for
% nonzeros we may pick. 

if ~with_replacement && (nsamps > nzrx)
    error('Cannot sample more than the total number of zeros');
end

% Save the requested number of zeros
nsamps_requested = nsamps;

% First, determine the number of samples we need to account for the fact
% that some samples will be nonzeros and so discarded.
ntmp1 = ceil(nsamps * nelx / nzrx);

% Error check
if ~with_replacement && (ntmp1 >= nelx)
    error('Need too many zero samples for this to work');
end

% Second, determine number of samples given that some will be duplicates,
% via coupon collector problem. This only matters if sampling with
% replacement.

if with_replacement
    ntmp2 = ntmp1;
else
    ntmp2 = ceil(nelx * log(1/(1-(ntmp1/nelx))));
end

% Finally, add a margin of safety by oversampling
nsamps = ceil(oversample * ntmp2); 

%% Generate the actual samples, removing duplicates, nonzeros, and excess

% Subscripts
tmpsubs = bsxfun(@(a,b)ceil(a.*b), rand(nsamps,d), sz);

if ~with_replacement
    tmpsubs = unique(tmpsubs,'rows','stable');
end

% Select out just the zeros
tmpidx = tt_sub2ind64(sz, tmpsubs);
iszero = ~builtin('_ismemberhelper',tmpidx,Xnzidx);
tmpsubs = tmpsubs(iszero,:);

% Trim back to desired number of samples
nsamps = min(size(tmpsubs,1), nsamps_requested);


%% Final return values
% Warn if too few entries
if (nsamps < nsamps_requested)
    warning('Unable to get the desired number of zero samples, %d (obtained) versus %d (requested)', nsamps, nsamps_requested);
end

% Allocate and fill subscript-value pairs.
subs = tmpsubs(1:nsamps,:);
