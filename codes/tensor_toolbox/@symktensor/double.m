function A = double(X)
%DOUBLE Convert a symktensor to a double array.
%
%   A = double(X) converts X to a standard multidimensional array.
%
%   See also SYMKTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isempty(X.lambda) % check for empty tensor
    A = [];
    return;
end

A = double(ktensor(X));
