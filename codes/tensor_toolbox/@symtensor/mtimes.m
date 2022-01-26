function C = mtimes(A,B)
%MTIMES tensor-scalar multiplication.
% 
%   C = MTIMES(A,B) is called for the syntax 'A * B' when A or B is a
%   symtensor and the other argument is a scalar.
% 
%   For symtensor-symtensor array multiplication, use TIMES or 'A .* B'.
% 
%   Examples
%   X = symtenrand([4,4,4])
%   W = 5 * X
%
%   See also SYMTENSOR, SYMTENSOR/TIMES
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%%
if isscalar(B)
    C = A;
    C.val = B * C.val;
    return;
end

if isscalar(A)
    C = B;
    C.val = A * C.val;
    return;
end

error('Mtimes only supports a symtensor times a scalar');





