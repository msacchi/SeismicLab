function ok = tt_subscheck(subs)
%TT_SUBSCHECK Checks for valid subscripts.
%
%   TT_SUBSCHECK(S) throws an error if S is not a valid subscript
%   array, which means that S is a matrix of real-valued, finite,
%   positive, integer subscripts.
%
%   X = TT_SUBSCHECK(S) returns true if S is a valid and false
%   otherwise.
%
%   See also TT_SIZECHECK, TT_VALSCHECK.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

%
% Includes improvements offered by Marcus Brubaker.

if isempty(subs)
    ok = true;
elseif ndims(subs) == 2 && isreal(subs) ...
        && all(isfinite(subs(:)) & subs(:) > 0) ...
        && isequal(subs,round(subs)) 
    ok = true;
else
    ok = false;
end

if ~ok && nargout == 0
    error('Subscripts must be a matrix of real positive integers');
end
