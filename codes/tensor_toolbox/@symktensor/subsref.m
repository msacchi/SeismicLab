function a = subsref(t,s)
%SUBSREF Subscripted reference for a symktensor.
%
%   Subscripted reference for a symtensor can be used to return the weight
%   vector, factor matrix, or expanded components of a symktensor.
%
%   Examples
%   S = symktensor(3, symtensor(@rand,4,3)); % <-- Declare a symktensor 
%   S.lambda % <-- Returns the weight vector.
%   S.X % <-- Returns the factor matrix.
%   S.M % <-- Returns the order (same as ndims(X)).
%   S(2,3,1) % <-- Calculates and returns that single element of X.
%
%   See also SYMKTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



switch s(1).type    
    case '.'
        switch s(1).subs
            case 'lambda'
                a = tt_subsubsref(t.lambda,s);
            case {'u','U','X','x'}
                a = tt_subsubsref(t.u,s);
            case {'m','M'}
                a = tt_subsubsref(t.m,s);
            otherwise
                error(['No such field: ', s(1).subs]);
        end
    case '()'
        if length(s.subs) == t.m %Needs to be polished.
            a = 0;
            for k = 1 : length(t.lambda)
                b = t.lambda(k);
                for i = 1 : length(s.subs)
                    b = b * t.u(s.subs{i},k);
                end
                a  = a + b;
            end
        else
            error('Incorrect index length');
        end
    otherwise
        error('Invalid subsref.');
end
