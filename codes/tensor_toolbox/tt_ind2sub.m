function subs = tt_ind2sub(siz,idx)
%TT_IND2SUB Multiple subscripts from linear indices.
%
%   SUBS = TT_IND2SUB(SIZ,INDS) returns that subscripts equivalent
%   to the linear indices in INDS for a tensor of size SIZ.
%
%   See also TT_SUB2IND, IND2SUB.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if isempty(idx)
    subs = [];
    return;
end

k = [1 cumprod(siz(1:end-1))];
n = length(siz);
idx=idx-1;
for i = n : -1 : 1 
    div=floor(idx/k(i));
    subs(:,i) = div+1;
    idx=idx-k(i)*div;
end
