function Z = ne(X,Y)
%NE Not equal (~=) for tensors.
%
%   See also TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



Z = tenfun(@ne,X,Y);
