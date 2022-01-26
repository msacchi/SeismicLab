function V = mttkrp(X,U,n)
%MTTKRP Matricized tensor times Khatri-Rao product for ktensor.
%
%   V = MTTKRP(X,U,n) efficiently calculates the matrix product of the
%   n-mode matricization of X with the Khatri-Rao product of all
%   entries in U, a cell array of matrices, except the nth.  How to
%   most efficiently do this computation depends on the type of tensor
%   involved.
%
%   See also TENSOR/MTTKRP, KTENSOR, KTENSOR/TTV
%
%   Examples
%   K = ktensor([2; 4], rand(3,2), rand(4,2), rand(5,2));
%   mttkrp(K, {rand(3,6), rand(4,6), rand(5,6)}, 3)
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



N = ndims(X);

if (n==1)
    R = size(U{2},2);
else
    R = size(U{1},2);
end

% Compute matrix of weights
W = repmat(X.lambda,1,R);
for i = [1:n-1,n+1:N]
  W = W .* (X.u{i}' * U{i});
end    

% Find each column of answer by multiplying columns of X.u{n} with weights 
V = X.u{n} * W;
