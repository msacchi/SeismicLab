function M = update(M,modes,data)
%UPDATE Update one or more modes of the ktensor with new data.
%
%   K = update(K,N,U) updates the ktensor K in mode N with the data in U
%   (in vector or matrix form). The value of N must be an integer between 1
%   and NDIMS(K). Further, NUMEL(U) must equal SZ(K,N) * NCOMPONENTS(K).
%
%   K = update(K,0,LAMBDA) updates the ktensor's lambda vector with the
%   data in LAMBDA, which must be a vector of length NDIMS(K).
%
%   K = update(K, MODES, DATA) updates the modes of K that are specified by
%   MODES with DATA. Here we assume that MODES is an ordered subset of 
%   {0,1,...,NDIMS(K)}. The vector DATA is a concatenation of vector data
%   to be used to replace each mode.  This mode is particularly useful when
%   working with an optimization method that understands only vectors of
%   unknowns.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

%% Error checking
if nargin < 3
    error('Update requires three arguments');
end

if ~isscalar(modes)
    if ~isequal(modes,sort(modes,'ascend'))
        error('Modes must be sorted');
    end
end

%%
loc = 1; % Location in data array
sz = size(M);
r = ncomponents(M);
for k = modes
    if k == 0
        endloc = loc + r - 1;
        if length(data) < endloc
            error('Data is too short');
        end
        M.lambda = data(loc:endloc);
        loc = endloc+1;        
    elseif k <= ndims(M)
        endloc = loc + sz(k)*r - 1;
        if length(data) < endloc
            error('Data is too short');
        end
        M.u{k} = reshape(data(loc:endloc),sz(k),r);
        loc = endloc+1;                        
    else
        error('Invalid mode: %d', k);
    end
    
end

%% Check that we used all the data
if loc ~= length(data)+1
    warning('Failed to consume all of the input data');
end

