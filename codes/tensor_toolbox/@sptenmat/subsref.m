function a = subsref(t,s)
%SUBSREF Subscripted reference for a sptenmat.
%
%   Examples
%   A.subs <-- returns the nonzero values as an array
%   A.vals <-- returns the corresponding 2D subscripts
%   A.tsize <-- returns the size original tensor
%   A.rdims <-- tensor dimensions that were mapped to rows
%   A.cdims <-- tensor dimensions that were mapped to columns 
%
%   See also SPTENMAT.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



switch s(1).type    
    case '.'
        switch s(1).subs
            case 'vals'
                a = tt_subsubsref(t.vals,s);
	    case 'tsize'
		a = t.tsize;
	    case 'rdims'
		a = t.rdims;
	    case 'cdims'
		a = t.cdims;
	    case 'subs' 
                a = tt_subsubsref(t.subs,s);
            otherwise
                error(['No such field: ', s.subs]);
        end
    otherwise
        error('Invalid subsref into tenmat.')
end
