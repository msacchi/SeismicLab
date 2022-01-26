function a = subsref(t,s)
%SUBSREF Subscripted reference for tenmat.
%
%   Examples
%   X(i,j) <-- returns the (i,j) entry in X
%   X.data <-- returns a 2D array of the data
%   X.tsize <-- returns the size original tensor
%   X.rdims <-- tensor dimensions that were mapped to rows
%   X.cdims <-- tensor dimensions that were mapped to columns 
%
%   See also TENMAT.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



switch s(1).type
    case '.'
        switch s(1).subs
            case 'data'
                a = tt_subsubsref(t.data,s);
            case 'tsize'
                a = tt_subsubsref(t.tsize,s);
            case {'rindices','rdims'}
                a = tt_subsubsref(t.rindices,s);
            case {'cindices','cdims'}
                a = tt_subsubsref(t.cindices,s);
            otherwise
                error(['No such field: ', s.subs]);
        end
    case '()'
        a = t.data(s.subs{:});
    otherwise
        error('Invalid subsref into tenmat.')
end
end
