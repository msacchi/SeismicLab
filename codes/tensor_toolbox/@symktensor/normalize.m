function X = normalize(X,N,normtype)
%NORMALIZE Normalizes the columns of the factor matrix.
%
%   NORMALIZE(X) normalizes the columns of the factor matrix using the
%   vector 2-norm, absorbing the excess weight into lambda. 
%
%   NORMALIZE(X,0) absorbs the weight into the factor matrix.
%   (All the lambda values are +/-1.)
%
%   NORMALIZE(X,[]) is equivalent to NORMALIZE(X). 
%
%   NORMALIZE(X,'sort') is the same as the above except it sorts the
%   components by lambda value, from greatest magnitude to least. 
%
%   NORMALIZE(X,V,1) normalizes using the vector one norm (sum(abs(x))
%   rather than the two norm (sqrt(sum(x.^2))), where V can be any of the
%   second arguments decribed above.
%
%   See also SYMKTENSOR, SYMKTENSOR/ARRANGE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%%
if ~exist('N','var')
    N = -1;
end

if isempty(N)
    N = -1;
end

if isequal(N,'sort')
    N = -2;
end

if N > 0
    error('Invalid second argument');
end

if ~exist('normtype','var')
    normtype = 2;
end

%% Ensure that matrix is normalized
for r = 1:length(X.lambda)
    tmp = norm(X.u(:,r),normtype);
    if (tmp > 0)
        X.u(:,r) = X.u(:,r) / tmp;
    end
    X.lambda(r) = X.lambda(r) * tmp.^(X.m);
    
    % Odd-ordered tensors should not have negative lambda values
    if (X.lambda(r) < 0) && (mod(X.m,2) == 1)
        X.u(:,r) = -X.u(:,r);
        X.lambda(r) = -X.lambda(r);
    end    
end

%% Absorb the weight into one factor, if requested
if (N == 0)
    d = nthroot(abs(X.lambda),X.m);
    X.u = bsxfun(@times,X.u',d)'; 
    X.lambda = sign(X.lambda) .* ones(size(X.lambda));
elseif (N == -2)
    if ncomponents(X) > 1
        [~,p] = sort(abs(X.lambda),'descend');
        X = arrange(X,p);
    end
end



