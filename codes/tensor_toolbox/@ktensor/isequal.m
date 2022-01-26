function [tf, tf_lambda, tf_U] = isequal(A,B)
%ISEQUAL True if each datum of two ktensor's are numerically equal.
%
%   TF = ISEQUAL(A,B) returns true if each factor matrix and the lambda
%   values are equal for A and B. Does not do any scaling or normalization
%   first.
%
%   [TF, TF_LAMBDA, TF_FACTORS] = ISEQUAL(A,B) returns also the result of
%   comparing the lambda vectors (TF_LAMBDA) and an array with the results
%   of comparing the factor matrices (TF_FACTORS).
%
%   See also KTENSOR, KTENSOR/NORMALIZE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


tf = false;
tf_lambda = false;
tf_U = false;

if ~isa(B,'ktensor')
    return;
end    

tf_lambda = isequal(A.lambda, B.lambda);
if ncomponents(A) == ncomponents(B)
    tf_U = cellfun(@isequal, A.u, B.u);
end
tf = tf_lambda & all(tf_U);

