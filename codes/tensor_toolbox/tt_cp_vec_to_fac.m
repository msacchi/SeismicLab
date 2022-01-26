function A = tt_cp_vec_to_fac(x,Z)
%TT_CP_VEC_TO_FAC Converts a vector to a cell array of factor matrices.
%
%   A = TT_CP_VEC_TO_FAC(X,Z) converts the vector X into a cell array
%   of factor matrices consistent with the size of the tensor Z.
%
%   See also TT_FAC_TO_VEC, TT_CP_FUN, TT_CP_OPT.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%% Set-up
P = length(x);
N = ndims(Z);
sz = size(Z);

%% Determine R
R = P / sum(sz);

%% Create A
A = cell(N,1);
for n = 1:N
    idx1 = sum(sz(1:n-1))*R + 1;
    idx2 = sum(sz(1:n))*R;
    A{n} = reshape(x(idx1:idx2),sz(n),R);
end
