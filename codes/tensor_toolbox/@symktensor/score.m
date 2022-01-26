function [best_score, A, flag, best_perm] = score(A,B,varargin)
%SCORE Checks if two symktensors match except for permutation.
%   
%   SCORE(A,B) returns the score of the match between A and B where
%   A is trying to be matched against B. It converts both to single-mode
%   ktensors and calls the ktensor SCORE function. 
%
%   See also SYMKTENSOR, KTENSOR/SCORE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%% Make sure A and B are symmetric ktensors
if ~isa(A,'symktensor') || ~isa(B,'symktensor')
    error('Both arguments must be symktensors');
end

A = normalize(A);
B = normalize(B);

AA = ktensor(A.lambda, A.u);
BB = ktensor(B.lambda, B.u);

[best_score, A, flag, best_perm] = score(AA,BB,varargin{:});
