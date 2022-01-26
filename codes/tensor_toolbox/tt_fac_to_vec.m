function x = tt_fac_to_vec(A)
%TT_FAC_TO_VEC Converts a set of factor matrices to a vector.
%
%   X = TT_FAC_TO_VEC(A) converts a cell array of factor matrices A to a
%   vector by vectorizing each matrix and stacking them.
%
%   See also TT_CP_VEC_TO_FAC, TT_CP_FUN, CP_OPT.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%% Set-up
N = length(A);

%% Get sizes
sz = zeros(N,1);
for n = 1:N
    sz(n) = size(A{n},1);
end
R = size(A{1},2);
P = sum(sz)*R;

%% Create x
x = zeros(P,1);
for n = 1:N
    idx1 = sum(sz(1:n-1))*R + 1;
    idx2 = sum(sz(1:n))*R;
    x(idx1:idx2) = reshape(A{n},sz(n)*R,1);
end
