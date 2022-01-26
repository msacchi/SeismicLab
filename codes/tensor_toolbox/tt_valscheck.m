function ok = tt_valscheck(vals)
%TT_VALSCHECK Checks for valid values.
%
%   TT_VALSCHECK(S) throws an error if S is not a valid values
%   array, which means that S is a column array.
%
%   X = TT_VALSCHECK(S) returns true if S is a valid and false otherwise.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if isempty(vals)
    ok = true;
elseif ndims(vals) == 2 && size(vals,2) == 1
    ok = true;
else
    ok = false;
end

if ~ok && nargout == 0
    error('Values must be a column array');
end
