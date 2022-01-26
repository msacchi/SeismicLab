function ok = tt_sizecheck(siz)
%TT_SIZECHECK Checks that the size is valid.
%
%   TT_SIZECHECK(S) throws an error if S is not a valid size array,
%   which means that it is a row vector with strictly postitive,
%   real-valued, finite integer values.
%
%   X = TT_SIZECHECK(S) returns true if S is a valid and false otherwise.
%
%   See also TT_SUBSCHECK.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if ndims(siz) == 2 && size(siz,1) == 1 ...
        && isreal(siz) ...
        && ~any(isnan(siz(:))) && ~any(isinf(siz(:))) ...
        && isequal(siz,round(siz)) && all(siz(:) > 0) 
    ok = true;
elseif isempty(siz)
    ok = true;    
else
    ok = false;
end

if ~ok && nargout == 0
    error('Size must be a row vector of real positive integers');
end
