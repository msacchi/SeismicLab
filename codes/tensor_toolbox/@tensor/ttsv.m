function y = ttsv(A,x,n,ver)
%TTSV Tensor times same vector in multiple modes.
%
%   Y = TTSV(A,X) multiples the tensor A by the vector X in all modes.
%
%   Y = TTSV(A,X,-1) multiplies the tensor A by the vector X in all modes
%   but the first. Returns the answer as a normal MATLAB array (not a
%   tensor). 
%
%   Y = TTSV(A,X,-2) multiplies the tensor A by the vector X in all modes
%   but the first two. Returns the answer as a normal MATLAB matrix (not a
%   tensor). 
%
%   Y = TTSV(A,X,-N) multiplies the tensor A by the vector X is all modes
%   but the first N.
%
%   See also TENSOR, TENSOR/TTV.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%% Process inputs (only two simple cases are supported)
if ~exist('n','var')
    n = 0;
elseif (n > 0)
    error('Invalid usage');
end

if ~exist('ver','var')
    ver = 2;
end


%% Calculate - old way
if ver == 1
    P = ndims(A);
    [X{1:P}] = deal(x);
    if (n == 0)
        y = ttv(A,X);
    elseif (n == -1) || (n == -2)
        y = double(ttv(A,X,-(1:-n)));
    else
        y = ttv(A,X,-(1:-n));
    end
    return
end

%% Calculate - new way
d = ndims(A); % Number of modes in tensor
sz = size(A,1); % Size of one mode (they're all the same)

dnew = -n; % Number of modes in result
drem = d - dnew; % Number of modes being multiplied out

%[X{1:drem}] = deal(x);
% Xkr = khatrirao(X);
% Ars = reshape(A.data, sz.^dnew, sz.^drem);
% y =  Ars * Xkr;

y = A.data;
for i = drem: -1 : 1
   yy = reshape(y,sz.^(dnew + i - 1),sz); 
   y = yy * x; 
end

% Convert to matrix if 2-way or convert back to a tensor if the result is
% 3-way or higher. Leave scalar or vector result alone.
if (dnew == 2)
    y = reshape(y, [sz sz]);
elseif dnew > 2
    y = tensor(y, sz * ones(dnew,1));
end


