function [I,C,W,Q] = indices(varargin)
%INDICES Compute unique indices of a symmetric tensor.
%
%   [I,C,W,Q] = INDICES(A) returns all unique indices for a
%   symmetric tensor. Each row of I is an index listed in increasing order.
%   Each row of C is the corresponding monomial representation, and W 
%   is the count of how many times that index appears in the symmetric 
%   tensor. Q is the number of rows of I, the number of unique indices.
%
%   See also SYMTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if nargin == 0
    error('INDICES requires at least one input argument');
elseif nargin == 2 % Specify m and n      
    m = varargin{1};
    n = varargin{2};    
elseif nargin == 1
    A = varargin{1};
    if ~(isa(A,'symktensor') || isa(A,'tensor') || isa(A,'symtensor'))
        error('First argument must be a scalar or a tensor-like class');
    end
    m = ndims(A);
    n = size(A,1);
else
    error('Wrong number of input arguments');
end

%% Determine size
sz = nchoosek(m+n-1,m);

%% Create I
% Following from function UpdateIndex (Figure 4) in
% G. Ballard, T. G. Kolda and T. Plantenga, Efficiently Computing Tensor
% Eigenvalues on a GPU, IPDPSW'11: Proceedings of the 2011 IEEE
% International Symposium on Parallel and Distributed Processing Workshops
% and PhD Forum, 12th IEEE International Workshop on Parallel and
% Distributed Scientific and Engineering Computing (PDSEC-11), Anchorage,
% Alaska (2011-05-16 to 2011-05-20), IEEE Computer Society, pp. 1340-1348,
% May 2011, doi:10.1109/IPDPS.2011.287       

I = zeros(sz,m);

for loc = 1:sz
    if loc == 1
        I(loc,:) = ones(1,m);
    else
        I(loc,:) = I(loc-1,:);
        j = m;
        while (I(loc,j) == n)
            j = j - 1;
        end
        I(loc,j:m) = I(loc,j)+1;
    end
end

if nargout==1    %Function can be called without monomials or weights
    return
end

%% Compute C from I
C = zeros(sz,n);
for i = 1:n
    C(:,i) = sum(I == i,2);
end

%% COMPUTE W (weights) from C
W = zeros(sz,1);
for i = 1:sz
    W(i) = multinomial(m,C(i,:));
end
Q=sz;
end
