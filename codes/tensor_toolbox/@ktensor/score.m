function [best_score, A, flag, best_perm] = score(A,B,varargin)
%SCORE Checks if two ktensors match except for permutation.
%   
%   SCORE(A,B) returns the score of the match between A and B where
%   A is trying to be matched against B.
%  
%   We define matching as follows. If A and B are single component ktensors 
%   that have been normalized so that their weights are lambda_a and
%   lambda_b, then the score is defined as
%   
%      score = penalty * (a1'*b1) * (a2'*b2) * ... * (aR'*bR),
%     
%   where the penalty is defined by the lambda values such that
%
%      penalty = 1 - abs(lambda_a - lambda_b) / max(lamdba_a, lambda_b).
%
%   The score of multi-components ktensors is a normalized sum of the
%   scores across the best permutation of the components of A. A can have
%   more components than B --- any extra components are ignored in terms of
%   the matching score.     
%
%   [SCORE,A] = SCORE(...) also returns A which has been normalized
%   and permuted to best match B. 
%
%   [SCORE,A,FLAG] = SCORE(...) also returns a boolean to indicate
%   a match according to a user-specified threshold.
%
%   [SCORE,A,FLAG,PERM] = SCORE(...) also returns the permutation
%   of the components of A that was used to best match B. 
%
%   SCORE(A,B,'param',value,...) takes the following parameters...
%
%      'lambda_penalty' - Boolean indicating whether or not to consider the
%      lambda values in the calculations. Default: true
%
%      'threshold' - Threshold specified in the formula above for
%      determining a match. Default: 0.99^N where N = ndims(A)
%
%      'greedy' - Boolean indicating whether or not to consider all
%      possible matchings (exponentially expensive) or just do a greedy
%      matching. Default: true
%
%   Examples
%   A = ktensor([2; 1; 2], rand(3,3), rand(4,3), rand(5,3));
%   B = ktensor([2; 4], ones(3,2), ones(4,2), ones(5,2));
%   score(A, B) %<--score(B,A) does not work: B has more components than A
%   score(A, B, 'greedy', false) %<--Check all permutations
%   score(A, B, 'lambda_penalty', false) %<--Without lambda penalty
%
%   This method is described in G. Tomasi and R. Bro, A Comparison of
%   Algorithms for Fitting the PARAFAC Model, Computational Statistics &
%   Data Analysis, Vol. 50, No. 7, pp. 1700-1734, April 2006,
%   doi:10.1016/j.csda.2004.11.013.
%  
%   See also KTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%  E. Acar & T. Kolda, 2010.

%% Make sure A and B are ktensors
if ~isa(A,'ktensor')
    A = ktensor(A);
end
if ~isa(B,'ktensor')
    B = ktensor(B);
end

%% Error checking
if ~isequal(size(A),size(B))
    error('Size mismatch');
end

%% Set-up
N = ndims(A);
RA = ncomponents(A);
RB = ncomponents(B);

%% We're matching components in A to B
if (RA < RB)
    error('Tensor A must have at least as many components as tensor B.');
end

%% Parse parameters
params = inputParser;
params.addParamValue('lambda_penalty', true, @islogical);
params.addParamValue('greedy', true, @islogical);
params.addParamValue('threshold', 0.99^N, @(x)(x<1));
params.parse(varargin{:});

%% Make sure columns of factor matrices in A and B are normalized
A = normalize(A);
B = normalize(B);

%% Compute all possible vector-vector congruences.

% Compute every pair for each mode
Cbig = tenzeros([RA,RB,N]);
for n = 1:N
    Cbig(:,:,n) = abs(A.u{n}' * B.u{n});
end

% Collapse across all modes using the product
C = double(collapse(Cbig,3,@prod));

%% Calculate penalty based on differences in the Lambda's
% Note that we are assuming the the lambda value are positive because the
% ktensor's were previously normalized.
if (params.Results.lambda_penalty)
    P = zeros(RA,RB);
    for ra = 1:RA
        la = A.lambda(ra);
        for rb = 1:RB
            lb = B.lambda(rb);
            if (la == 0) && (lb == 0)
                % if both lambda values are zero (0), they match
                P(ra, rb) = 1;
            else
                P(ra,rb) = 1 - (abs(la-lb) / max(abs(la),abs(lb)));
            end
        end
    end
    C = P.*C;
end

%% Option to do greedy matching
if (params.Results.greedy)
    
    best_perm = zeros(1,RA);
    best_score = 0;
    for r = 1:RB
        [~,idx] = max(C(:));
        [i,j] = ind2sub([RA RB], idx);
        best_score = best_score + C(i,j);
        C(i,:) = -10;
        C(:,j) = -10;
        best_perm(j) = i;
    end
    best_score = best_score / RB;
    flag = 1;
    
    % Rearrange the components of A according to the best matching
    foo = 1:RA;
    tf = ismember(foo,best_perm);
    best_perm(RB+1:RA) = foo(~tf);
    A = arrange(A, best_perm);
    return;
end

%% Compute all possible matchings
% Creates a matrix P where each row is a possible matching of components in
% A to components of B. We assume A has at least as many components as B.
idx = nchoosek(1:RA,RB);
M = [];
for i = 1:size(idx,1)
    M = [M; perms(idx(i,:))]; %#ok<AGROW>
end

%% Calculate the congruences for each matching
scores = zeros(size(M));
for i = 1:size(M,1)
    for r = 1:RB
        scores(i,r) = C(M(i,r),r);
    end
end

%% Figure out the best matching based on sum's across the components
score = sum(scores,2)/RB;
[best_score, max_score_id] = max(score);
if min(scores(max_score_id,:)) >= params.Results.threshold
    flag = 1;
else
    flag = 0;
end
best_match = M(max_score_id,:);
best_perm = [best_match setdiff(1:RA, best_match)];

%% Rearrange the components of A according to the best matching
A = arrange(A, best_perm);


