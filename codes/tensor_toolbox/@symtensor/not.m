function B = not(A)
%NOT Logical NOT (~) for symmetric tensors.
%
%   See also SYMTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


B = A;
B.val = not(A.val);
