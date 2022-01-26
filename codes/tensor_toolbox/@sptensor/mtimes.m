function C = mtimes(A,B)
%MTIMES sptensor-scalar multiplication.
% 
%   C = MTIMES(A,B) is called for the syntax 'A * B' when A or B is a
%   sparse tensor and the other argument is a scalar.
% 
%   For tensor-matrix multiplication, use TTM.
%   For tensor-tensor multiplication, use TTT.
%   For tensor-tensor array multiplication, use TIMES or 'A .* B'.
%
%   See also SPTENSOR, SPTENSOR/TTM, SPTENSOR/TTT, SPTENSOR/TIMES
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isscalar(B)    
    C = sptensor(A.subs, A.vals * B, size(A));
    return;
end

if isscalar(A)
    C = sptensor(B.subs, B.vals * A, size(B));
    return;
end

error('MTIMES only supports the scalar case for sparse tensors');
