function X = full(T)
%FULL Convert a sumtensor to a (dense) tensor.
%
%   X = FULL(T) converts sumtensor T to (dense) tensor X. This may be an
%   expensive operation for large-scale tensors.
%
%   See also SUMTENSOR, TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if isempty(T.part)
    X = tensor;
    return;
end

X = full(T.part{1});
for i = 2:length(T.part)
    X = X + full(T.part{i});
end

return;
