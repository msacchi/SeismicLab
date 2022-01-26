function C = mtimes(A,B)
%MTIMES Implement A*B (scalar multiply) for ktensor.
%
%   C = mtimes(A,B) computes A * B where A is a Kruskal tensor and B is
%   a scalar (or vice versa). The result C is the same size as A.
%
%   See also KTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



% Note: We can do scalar times a tensor, but anything more complex is
% an error.

if isa(B,'numeric') && isequal(size(B),[1 1])
    C = ktensor(B * A.lambda, A.u);
elseif isa(A,'numeric') && isequal(size(A),[1 1])
    C = ktensor(A * B.lambda, B.u);
else
    error('Use mtimes(full(A),full(B)).');
end
