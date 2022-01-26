function C = mtimes(A,B)
%MTIMES tensor-scalar multiplication.
% 
%   C = MTIMES(A,B) is called for the syntax 'A * B' when A or B is a
%   tensor and the other argument is a scalar.
% 
%   For tensor-matrix multiplication, use TTM.
%   For tensor-tensor multiplication, use TTT.
%   For tensor-tensor array multiplication, use TIMES or 'A .* B'.
% 
%   Examples
%   X = tenrand([3,4,2])
%   W = 5 * X
%
%   See also TENSOR, TENSOR/TTM, TENSOR/TTT, TENSOR/TIMES
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%%
if isscalar(B)
    C = A;
    C.data = B * C.data;
    return;
end

if isscalar(A)
    C = B;
    C.data = A * C.data;
    return;
end

error('Mtimes only supports a tensor times a scalar');





