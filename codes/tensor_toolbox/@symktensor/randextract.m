function [Mnew,info] = randextract(M,s)
%RANDEXTRACT Create new symktensor with random subset of components.
%
%   MS = RANDEXTRACT(M,S) created a new matrix MS that comprises S
%   random columns of M. Note that columns from M may be repeated in MS.
%
%   MS = RANDEXTRACT(M,C) creates a new matrix MS by extracting the
%   columns specified in the vector C. 
%
%   MS = RANDEXTRACT(M,0) just returns M, i.e., no sampling. 
%
%   See also CP_ISYM.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%% Extract sizes
p = size(M.u,2);

%% Sample indices (or use existing sample)
if isequal(s,0)
    Mnew = M;
    return;
elseif isscalar(s)
    nidx = randi(p,1,s);
else
    nidx = s;
    s = length(nidx);
end

%% Extract from U
Unew = M.u(:,nidx);

%% Update lambda
Lnew = (p/s) * M.lambda(nidx);

%% Create new M
d = ndims(M);
Mnew = symktensor(Lnew,Unew,d);

%% Save info
if nargout == 2
    info.nidx = nidx;
end
