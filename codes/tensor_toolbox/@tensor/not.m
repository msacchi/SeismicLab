function B = not(A)
%NOT Logical NOT (~) for tensors.
%
%   See also TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



B = tensor(not(A.data), size(A));
