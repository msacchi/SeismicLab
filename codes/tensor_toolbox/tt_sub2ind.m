function idx = tt_sub2ind(siz,subs)
%TT_SUB2IND Converts multidimensional subscripts to linear indices.
%
%   INDS = TT_SUB2IND(SIZ,SUBS) returns the linear indices
%   equivalent to the subscripts in the array SUBS for a tensor of
%   size SIZ.  
%
%   See also TT_IND2SUB, SUB2IND.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if isempty(subs)
    idx = [];
    return;
end

mult = [1 cumprod(siz(1:end-1))];
idx = (subs - 1) * mult' + 1;

