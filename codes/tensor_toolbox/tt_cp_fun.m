function [f,g] = tt_cp_fun(x,Z,Znormsqr)
%TT_CP_FUN Calculate function and gradient for CP fit function.
%
%   [F,G] = TT_CP_FUN(X,Z) where X is a vector containing the entries of the
%   components of the model and Z is the tensor to be fit.
%
%   See also TT_CP_VEC_TO_FAC, TT_FAC_TO_VEC, TT_CP_FG, CP_OPT
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%% Convert x to a cell array of matrices
A = tt_cp_vec_to_fac(x,Z);

%% Call cp_fit and cp_gradient using cp_fg
[f,G] = tt_cp_fg(Z,A,Znormsqr);

%% Convert a cell array to a vector
g = tt_fac_to_vec(G);


