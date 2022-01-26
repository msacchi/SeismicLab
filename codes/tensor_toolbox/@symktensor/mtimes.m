function C = mtimes(A,B)
%MTIMES Implement A*B (scalar multiply) for symktensor.
%
%   C = mtimes(A,B) computes A * B where A is a symktensor and B is
%   a scalar (or vice versa). The result C is the same size as A.
%
%   See also SYMKTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


% Note: We can do scalar times a tensor, but anything more complex is
% an error.

if isa(B,'numeric') && isequal(size(B),[1 1])
    C = symktensor(B * A.lambda, A.u, A.m);
elseif isa(A,'numeric') && isequal(size(A),[1 1])
    C = symktensor(A * B.lambda, B.u, B.m);
else
    error('Use mtimes(full(A),full(B)).');
end
