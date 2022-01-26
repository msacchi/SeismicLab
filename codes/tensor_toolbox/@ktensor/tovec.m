function x = tovec(K,lambdaflag)
%TOVEC Convert Ktensor to vector.
%
%   V = TOVEC(K) converts the Ktensor to a column vector of length
%   sum(size(K)+1)*ncomponents(K). The format is 
%      [ K.lambda; K.U{1}(:); K.U{2}(:); ... ]
%
%   V = TOVEC(K,false) ignores lambda in the conversion, so the vector V is
%   of length P where = sum(size(K))*ncomponents(K).
%
%   Examples
%   K = ktensor([3; 2], rand(4,2), rand(5,2), rand(3,2));
%   V = tovec(K);
%   Kcopy = ktensor(V, size(K), ndims(K), true);
%   norm(K-Kcopy) %<- Zero (or close to it)
%
%   K = ktensor({rand(4,2), rand(5,2), rand(3,2)});
%   V = tovec(K);
%   Kcopy = ktensor(V, size(K), ndims(K), false);
%   norm(K-Kcopy) %<- Zero (or close to it)
%
%   See also KTENSOR, KTENSOR/SIZE, KTENSOR/NCOMPONENTS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

if ~exist('lambdaflag','var')
    lambdaflag = true;
end

xcell = cellfun(@(x) x(:), K.u, 'UniformOutput', false);
x = cell2mat(xcell);

if lambdaflag
    x = [K.lambda; x];
end
