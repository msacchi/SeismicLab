function C = times(A,B)
%TIMES Array multiplication for sparse tensors.
%
%   TIMES(A,B) is called for the syntax 'A .* B' when A or B is a 
%   sparse tensor. A and B must have the same size, unless one is a scalar.
%   A scalar can be multiplied by a sparse tensor of any size. 
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%% Observations for sparse matrix case.
% The result of a .* 5 is sparse.
% The result of a .* 0 is sparse.
% The result of a .* full(a) is sparse.

%%
if isscalar(B)
    C = sptensor(A.subs, A.vals * B, size(A));
    return;
end

if isscalar(A)
    C = sptensor(B.subs, B.vals * A, size(B));
    return;
end

if ~isequal(size(A),size(B))
    error('Must be two tensors of the same size');
end

switch class(B)
    case {'sptensor'}
        [csubs,ia,ib] = intersect(A.subs,B.subs,'rows');
        cvals = A.vals(ia) .* B.vals(ib);
        C = sptensor(csubs, cvals, size(A));
        return;
    case {'tensor'}
        csubs = A.subs;
        cvals = A.vals .* B(csubs); 
        C = sptensor(csubs, cvals, size(A));
        return;       
    case {'ktensor'}    
        csubs = A.subs;
        cvals = zeros(size(A.vals));       
        R = numel(B.lambda);
        N = ndims(A);
        for r = 1:R
            tvals = B.lambda(r) * A.vals;
            for n = 1:N
                v = B{n}(:,r);
                tvals = tvals .* v(csubs(:,n));
            end
            cvals = cvals + tvals;
        end
        C = sptensor(csubs, cvals, size(A));
        return;       
    otherwise
        error('Invalid second argument for sptensor/times');
end
