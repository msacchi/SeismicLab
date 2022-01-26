function x = tovec(S,nolambda)
%TOVEC Convert symktensor to vector representation.
%
%   V = TOVEC(S) converts a symktensor to a vector. It stacks the LAMBDA
%   vector on top of a vectorized version of the matrix X.
%
%   V = TOVEC(S,TRUE) just returns a vectorized version of the matrix
%   X. It requires LAMBDA=1.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if exist('nolambda','var') && nolambda 
    if any(S.lambda ~= 1)
        error('Not all lambda values are 1.')
    end
    x = reshape(S.u,[],1);
else
    x = [S.lambda; reshape(S.u,[],1)];
end
