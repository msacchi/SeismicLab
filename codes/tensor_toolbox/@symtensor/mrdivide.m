function Z = mrdivide(X,Y)
%MRDIVIDE Slash right division for symmetric tensors.
%
%   MRDIVIDE(A,B) is called for the syntax 'A / B' when A is a symmetric 
%   tensor and B is a scalar. 
%
%   Example
%   X = tenrand([4 4 4]);
%   X / 3
%
%   See also SYMTENSOR, SYMTENSOR/RDIVIDE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isscalar(Y)
    Z = tenfun(@rdivide,X,Y);
    return;
end

error('MRDIVIDE only supports the scalar case for symmetric tensors');

