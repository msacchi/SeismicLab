function A = double(X)
%DOUBLE Convert a ktensor to a double array.
%
%   A = double(X) converts X to a standard multidimensional array.
%
%   See also KTENSOR, KTENSOR/FULL.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isempty(X.lambda) % check for empty tensor
    A = [];
    return;
end

sz = [size(X) 1];
A = X.lambda' * khatrirao(X.u,'r')';
A = reshape(A,sz);
