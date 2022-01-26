function nrm = norm(A)
%NORM Frobenius norm of a symktensor.
%
%   NORM(T) returns the Frobenius norm of a symktensor.
%
%   See also SYMKTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



% Retrieve the factors of A
U = A.u;
UtU = A.u'*A.u;

% Compute the matrix of correlation coefficients
coefMatrix = A.lambda * A.lambda';
coefMatrix = coefMatrix .* ((UtU).^(A.m));

nrm = sqrt(abs(sum(coefMatrix(:))));

return;
