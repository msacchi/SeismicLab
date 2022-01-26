function [tf, tf_lambda, tf_U] = isequal(A,B)
%ISEQUAL True if each component of two symktensors is numerically equal.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


tf = false;
tf_lambda = false;
tf_U = false;

if ~isa(B,'symktensor')
    return;
end    

tf_lambda = isequal(A.lambda, B.lambda);
tf_U = isequal(A.u, B.u);
tf = tf_lambda & tf_U;

