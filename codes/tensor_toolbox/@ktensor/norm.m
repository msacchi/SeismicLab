function nrm = norm(A)
%NORM Frobenius norm of a ktensor.
%
%   NORM(T) returns the Frobenius norm of a ktensor.
%
%   See also KTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



% Retrieve the factors of A
U = A.u;

% Compute the matrix of correlation coefficients
coefMatrix = A.lambda * A.lambda';
for i = 1:ndims(A)
  coefMatrix = coefMatrix .* (U{i}'*U{i});
end

nrm = sqrt(abs(sum(coefMatrix(:))));

return;
