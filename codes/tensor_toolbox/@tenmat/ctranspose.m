function a = ctranspose(a)
%CTRANSPOSE Complex conjugate transpose for tenmat.
%
%   C = CTRANSPOSE(A) swaps the row and column indices of A. 
%
%   See also TENMAT.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



tmp = a.rindices;
a.rindices = a.cindices;
a.cindices = tmp;
a.data = ctranspose(a.data);
