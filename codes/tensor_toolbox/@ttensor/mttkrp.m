function V = mttkrp(X,U,n)
%MTTKRP Matricized tensor times Khatri-Rao product for ttensor.
%
%   V = MTTKRP(X,U,n) efficiently calculates the matrix product of the
%   n-mode matricization of X with the Khatri-Rao product of all
%   entries in U, a cell array of matrices, except the nth.  How to
%   most efficiently do this computation depends on the type of tensor
%   involved.
%
%   Examples
%   T = ttensor(tensor(rand(2,3,4)), {rand(5,2), rand(2,3), rand(4,4)});
%   mttkrp(T , {rand(5,6), rand(2,6), rand(4,6)}, 3)
%
%   See also TENSOR/MTTKRP, TTENSOR, TTENSOR/TTV
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


N = ndims(X);

if (n==1)
    R = size(U{2},2);
else
    R = size(U{1},2);
end

% Compute cell array of weights to multiply into core
W = cell(N,1);
for i = [1:n-1,n+1:N]
  W{i} = (X.u{i}' * U{i});
end    
Y = mttkrp(X.core,W,n);

% Find each column of answer by multiplying columns of X.u{n} with weights 
V = X.u{n} * Y;
