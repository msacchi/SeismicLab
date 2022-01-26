function V = mttkrp(X,U,n,vers)
%MTTKRP Matricized tensor times Khatri-Rao product for tensor.
%
%   V = MTTKRP(X,U,N) efficiently calculates the matrix product of the
%   n-mode matricization of X with the Khatri-Rao product of all
%   entries in U, a cell array of matrices, except the Nth.  How to
%   most efficiently do this computation depends on the type of tensor
%   involved.
%
%   V = MTTKRP(X,K,N) instead uses the Khatri-Rao product formed by the
%   matrices and lambda vector stored in the ktensor K. As with the cell
%   array, it ignores the Nth factor matrix. The lambda vector is absorbed
%   into one of the factor matrices.
%
%   NOTE: Updated to use BSXFUN per work of Phan Anh Huy. See Anh Huy Phan,
%   Petr Tichavský, Andrzej Cichocki, On Fast Computation of Gradients for
%   CANDECOMP/PARAFAC Algorithms, arXiv:1204.1586, 2012.
%
%   Examples
%   mttkrp(tensor(rand(3,3,3)), {rand(3,3), rand(3,3), rand(3,3)}, 2)
%   mttkrp(tensor(rand(2,4,5)), {rand(2,6), rand(4,6), rand(5,6)}, 3)
%
%   See also TENSOR, TENMAT, KHATRIRAO
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


% Multiple versions supported...
if ~exist('vers','var')
    vers = 1;
end

N = ndims(X);
if (N < 2)
    error('MTTKRP is invalid for tensors with fewer than 2 dimensions');
end

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

if ~iscell(U)
    error('Second argument should be a cell array or a ktensor');
end

if (length(U) ~= N)
    error('Cell array is the wrong length');
end

if n == 1
    R = size(U{2},2);
else
    R = size(U{1},2);
end
   
for i = 1:N
   if i == n, continue; end
   if (size(U{i},1) ~= size(X,i)) || (size(U{i},2) ~= R)
       error('Entry %d of cell array is wrong size', i);
   end
end

%% Computation

if vers == 0 % Old version of the code
    Xn = permute(X,[n 1:n-1,n+1:N]);
    Xn = reshape(Xn.data, size(X,n), []);
    Z = khatrirao(U{[1:n-1,n+1:N]},'r');
    V = Xn*Z;
    return;
end

szl = prod(size(X,1:n-1)); %#ok<*PSIZE>
szr = prod(size(X,n+1:N));
szn = size(X,n);

if n == 1
    Ur = khatrirao(U{2:N},'r');
    Y = reshape(X.data,szn,szr);
    V =  Y * Ur;
elseif n == N
    Ul = khatrirao(U{1:N-1},'r');
    Y = reshape(X.data,szl,szn);
    V = Y' * Ul;
else
    Ul = khatrirao(U{n+1:N},'r');
    Ur = reshape(khatrirao(U{1:n-1},'r'), szl, 1, R);
    Y = reshape(X.data,[],szr);
    Y = Y * Ul;
    Y = reshape(Y,szl,szn,R);
    if vers == 2
        V = bsxfun(@times,Ur,Y);
        V = reshape(sum(V,1),szn,R);
    else % default (vers == 1)
        V = zeros(szn,R);
        for r =1:R
            V(:,r) = Y(:,:,r)'*Ur(:,:,r);
        end
    end
end

