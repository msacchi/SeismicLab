function a = subsref(t,s)
%SUBSREF Subscripted reference for a ttensor.
%
%   Subscripted reference is used to query the components of a ttensor.
%
%   Examples
%   core = tensor(rand(2,2,2));
%   X = ttensor(core, rand(4,2), rand(5,2), rand(3,2));
%   X.core %<-- returns core array
%   X.U %<-- returns a cell array of three matrices
%   X.U{1} %<-- returns the matrix corresponding to the first mode.
%
%   See also TTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


switch s(1).type    
    case '.'
        switch s(1).subs
            case {'core','lambda'}
                a = tt_subsubsref(t.core,s);
            case {'U','u'}
                a = tt_subsubsref(t.u,s);
            otherwise
                error(['No such field: ', s.subs]);
        end
    case '()'
	error('Subsref with () not supported for ttensor.');
    case '{}'
	new_s(1).type = '.';
	new_s(1).subs = 'u';
	new_s(2:length(s)+1) = s;
	a = subsref(t, new_s);
     otherwise
        error('Invalid subsref.');
end
