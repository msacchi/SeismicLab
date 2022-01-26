function V = mttkrps(X,U)
%MTTKRPS Sequence of MTTKRP calculations for a tensor.
%
%   V = MTTKRPS(X,U) computes a cell array V such that 
%   V{k} = mttkrp(X, U, k) for k = 1,...,ndims(X). 
%
%   See also MTTKRP.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

% Written by J. Duersch, 2018.
%
% Copyright (2018) Sandia Corporation. Under the terms of Contract

% Obtain dimensions and optimal splitting.
sz = size(X);
d = length(sz);
s = min_split(sz);

% Output sequence V{k} = mttkrp(X,U,k)
V = cell(d,1);

% KRP over modes s+1:d.
K = khatrirao(U{s+1:d},'r');
% Partial MTTKRP with remaining modes <m1, m2, ... , ms, C>
W = reshape(X.data,[],size(K,1)) * K;

for k=1:s-1
    % Loop entry invariant: W has modes <mk, ... , ms, C>
    V{k} = mttv_mid(W, U(k+1:s));
    % Satisfy invariant.
    W = mttv_left(W, U{k});
end

% Exit state: W has modes <ms, C>.
V{s} = W;

% KRP over modes 1:s.
K = khatrirao(U{1:s},'r');
% Partial MTTKRP with remaining modes <m{s+1}, m{s+2}, ... , md, C>
W = reshape(X.data,size(K,1),[])' * K;

for k=s+1:d-1
    % Loop entry invariant: W has modes <mk, ... , md, C>
    V{k} = mttv_mid(W, U(k+1:d));
    % Satisfy invariant.
    W = mttv_left(W, U{k});
end

% Exit state: W has modes <md, C>.
V{d} = W;
end




function W_out = mttv_left(W_in, U1)
% W_out = mttv_left(W_in, U_left)
% Contract leading mode in partial MTTKRP W_in using the matching factor
% matrix U1. The leading mode is defined as the mode for which consecutive
% increases in the corresponding index address elements at consecutive
% increases in the memory offset.
% 
% W_in has modes in natural descending order: <m1, m2, ... , mN, C>.
%    Mode m1 is either the first mode or an intermediate mode of the
%    original tensor. Mode m2 through mN are subsequence original modes.
%    The last mode C is the component mode (indexed over rank-1 components
%    1:r) corresponding to columns in factor matrices.
% U1 is the corresponding factor matrix with modes <m1, C>.
% W_out has modes: <m2, ... , mN, C>

r = size(U1,2);
W_in = reshape(W_in, size(U1,1), [], r);
W_out = zeros(size(W_in,2), r);
for j=1:r
    W_out(:,j) = W_in(:,:,j)' * U1(:,j);
end
end




function V = mttv_mid(W_in, U_mid)
% V = mttv_mid(W_in, U_mid)
% Contract all intermediate modes in partial MTTKRP W_in using the matching
% cell array U_mid.
% 
% W_in has modes in natural descending order: <m1, m2, ... , mN, C>.
%    Mode m1 is either the first mode or an intermediate mode of the
%    original tensor. Mode m2 through mN are subsequence original modes.
%    The last mode C is the component mode (indexed over rank-1 components
%    1:r) corresponding to columns in factor matrices.
% U_mid is the corresponding cell array of factor matrices. That is,
%    U_mid{1} has modes <m2, C>, U_mid{2} has modes <m3, C>, etc. The cell
%    array must exactly match all intermediate uncontracted modes.
% V is the final MTTKRP with modes: <m1, C>.
if isempty(U_mid)
    V = W_in;
else
    K = khatrirao(U_mid,'r');
    r = size(K,2);
    W_in = reshape(W_in, [], size(K,1), r);
    V = zeros(size(W_in,1), r);
    for j=1:r
        V(:,j) = W_in(:,:,j)*K(:,j);
    end
end
end




function [s_min]=min_split(sz)
% [s_min]=min_split(sz)
% Scan for optimal splitting with minimal memory footprint.
%
% sz gives sizes of each dimension in the original tensor in natural
%    descending order.
% s_min gives optimal splitting to minimize partial MTTKRP memory
%    footprint. Modes 1:s_min will contract in left-partial computation and
%    modes s_min+1:d will contract in right-partial.

m_left=sz(1);
m_right=prod(sz(2:end));
s_min=1;

% Minimize: m_left + m_right.
for s=2:length(sz)-1
    % Peel mode s off right and test placement.
    m_right = m_right/sz(s);
    if (m_left < m_right)
        % The sum is reduced by placing mode s on the left.
        s_min = s;
        m_left = m_left*sz(s);
    else
        % The sum would be reduced by placing mode s back on the right.
        % There is no further benefit to collecting modes on the left.
        break;
    end
end
end


