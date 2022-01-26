function C = xor(A,B)
%XOR Logical XOR for sptensors.
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%% Observations for sparse matrix case.
% The result of xor(a,5) is dense!
% The result of xor(a,0) is dense!
% The result of xor(a,full(a)) is dense!
% The result of xor(a,b) is sparse.

%% Case 1: One argument is a scalar
if isscalar(B) || isa(B,'tensor')
    C = xor(full(A),B);
    return;
end
if isscalar(A)
    C = xor(A,full(B));
    return;
end


%% Case 2: Both x and y are tensors of some sort
if ~isequal(size(A),size(B))
    error('Must be tensors of the same size');
end

if isa(A,'sptensor') && isa(B,'sptensor')
    C = sptensor([A.subs; B.subs], 1, size(A), @(x) length(x) == 1);
    return;
end

if isa(B,'tensor')
    Bsubs = find(B ~= 0);
    C = sptensor([A.subs; Bsubs], 1, size(A), @(x) length(x) == 1);
    return;    
end

%% Otherwise
error('The arguments must be two sptensors or an sptensor and a scalar.');
