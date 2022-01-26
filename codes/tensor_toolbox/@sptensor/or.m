function C = or(A,B)
%OR Logical OR (|) for sptensors.
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%% Observations for sparse matrix case.
% The result of a | 5 is dense!
% The result of a | 0 is dense!
% The result of a | full(a) is dense!
% The result of a | a is sparse.

%% Case 1: One argument is a scalar
if isscalar(B) || isa(B,'tensor')
    C = full(A) | B;
    return;
end
if isscalar(A)
    C = A | full(B);
    return;
end

%% Case 2: Both A and B are sparse tensors
if ~isequal(size(A),size(B))
    error('Must be tensors of the same size');
end

if isa(A,'sptensor') && isa(B,'sptensor')
    C = sptensor([A.subs; B.subs], 1, size(A), @(x) length(x) >= 1);
    return;
end

%% Otherwise
error('The arguments must be two sptensors or an sptensor and a scalar.');
