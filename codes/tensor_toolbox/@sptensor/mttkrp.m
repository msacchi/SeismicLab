function V = mttkrp(X,U,n)
%MTTKRP Matricized tensor times Khatri-Rao product for sparse tensor.
%
%   V = MTTKRP(X,U,n) efficiently calculates the matrix product of the
%   n-mode matricization of X with the Khatri-Rao product of all
%   entries in U, a cell array of matrices, except the nth.  How to
%   most efficiently do this computation depends on the type of tensor
%   involved.
%
%   V = MTTKRP(X,K,N) instead uses the Khatri-Rao product formed by the
%   matrices and lambda vector stored in the ktensor K. As with the cell
%   array, it ignores the Nth factor matrix. The lambda vector is absorbed
%   into one of the factor matrices.
%
%   Examples
%   S = sptensor([3 3 3; 1 3 3; 1 2 1], 4, [3, 4, 3]); %<-Declare sptensor
%   mttkrp(S, {rand(3,3), rand(3,3), rand(3,3)}, 2)
%
%   See also TENSOR/MTTKRP, SPTENSOR/TTV, SPTENSOR
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



% In the sparse case, it is most efficient to do a series of TTV operations
% rather than forming the Khatri-Rao product.

N = ndims(X);

if isa(U,'ktensor')
    % Absorb lambda into one of the factors, but not the one that's skipped
    if n == 1
        U = redistribute(U,2);
    else
        U = redistribute(U,1);
    end    
    % Extract the factor matrices
    U = U.u;
end

if (length(U) ~= N)
    error('Cell array is the wrong length');
end

if ~iscell(U)
    error('Second argument should be a cell array or a ktensor');
end

if (n == 1)
    R = size(U{2},2);
else
    R = size(U{1},2);
end

V = zeros(size(X,n),R);
for r = 1:R
    % Set up cell array with appropriate vectors for ttv multiplication
    Z = cell(N,1);
    for i = [1:n-1,n+1:N]
        Z{i} = U{i}(:,r);
    end
    % Perform ttv multiplication
    V(:,r) = double(ttv(X, Z, -n));
end
