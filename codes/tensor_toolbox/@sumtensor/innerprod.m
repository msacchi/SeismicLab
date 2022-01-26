function res = innerprod(X,Y)
%INNERPROD Efficient inner product with a sumtensor.
%
%   R = INNERPROD(X,Y) efficiently computes the inner product between
%   two tensors X and Y, where X is a sumtensor.
%
%   See also TENSOR/INNERPROD, SPTENSOR/INNERPROD, TTENSOR/INNERPROD, 
%   KTENSOR/INNERPROD
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


% X is a sumtensor

res = innerprod(X.part{1},Y);
for i = 2:length(X.part)
    res = res + innerprod(X.part{i},Y); 
end



