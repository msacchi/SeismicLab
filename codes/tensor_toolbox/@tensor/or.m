function Z = or(X,Y)
%OR Logical OR (|) for tensors.
%
%   See also TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



Z = tenfun(@or,X,Y);
