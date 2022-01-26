function [tf, tf_core, tf_U] = isequal(A,B)
%ISEQUAL True if the part of two ttensor's are numerically equal.
%
%   TF = ISEQUAL(A,B) returns true if each factor matrix and the core
%   are equal for A and B.  
%
%   [TF, TF_CORE, TF_FACTORS] = ISEQUAL(A,B) returns also the result of
%   comparing the core (TF_CORE) and an array with the results of comparing
%   the factor matrices (TF_FACTORS).
%
%    See also TTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


tf = false;
tf_core = false;
tf_U = false;

if ~isa(B,'ttensor')
    return;
end

if ndims(A) ~= ndims(B)
    return;
end

tf_core = isequal(A.core, B.core);
tf_U = cellfun(@isequal, A.u, B.u);
tf = tf_core & all(tf_U);

