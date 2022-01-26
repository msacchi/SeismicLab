function Y = divide(X,K,epsilon)
%DIVIDE Divide an SPTENSOR by a nonnegative KTENSOR.
%
%   Y = DIVIDE(X,K,EPSILON) divides the sparse tensor X by the 
%   nonnegative ktensor K.  Avoids divide-by-zero errors by dividing 
%   by MIN(EPSILON,K-VALUE) at each nonzero of X.
%
%   See also SPTENSOR, KTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


% Assumes K is a nonnegative ktensor

Y = X;

subs = Y.subs;
vals = zeros(size(Y.vals));
R = numel(K.lambda);
N = ndims(Y);
for r = 1:R
    tvals = ones(size(vals)) * K.lambda(r);
    for n = 1:N
        v = K{n}(:,r);
        tvals = tvals .* v(subs(:,n));
    end
    vals = vals + tvals;
end
Y.vals = Y.vals ./ max(epsilon, vals);

return;
