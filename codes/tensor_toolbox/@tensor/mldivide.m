function Z = mldivide(X,Y)
%MLDIVIDE Slash left division for tensors.
%
%   MLDIVIDE(A,B) is called for the syntax 'A \ B' when A is a scalar and B
%   is a tensor.  
%
%   Example
%   X = tenrand([4 3 2],5);
%   3 \ X
%
%   See also TENSOR, TENSOR/LDIVIDE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isscalar(X)
    Z = tenfun(@ldivide,X,Y);
    return;
end

error('MLDIVIDE only supports the scalar case for tensors');

