function t = subsasgn(t,s,b)
%SUBSASGN Subscripted assignment for symktensor.
%
%   Subscript assignment can be used to change the order of a symmetric
%   tensor decomposition, the weight vector, or the factor matrix.
%
%   Examples
%   S = symktensor(2, symtensor('rand', 4,3)); %<--Declare a symtensor
%   S.lambda = [2; 1]; %<-- Change the weight vector
%   S.X = rand(3,2); %<-- Change the factor matrix
%   S.U = rand(3,2); %<-- Same as above
%   S.m = 5; %<-- Change the order of the decomposition.
%
%   See also SYMKTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



switch s(1).type
    case '.'
        switch s(1).subs
            case {'m','M'}
                if isscalar(b)
                    t.m = b;
                else
                    error('M must be a scalar');
                end
            case 'lambda'
                if length(b) == length(t.lambda)
                    t.lambda = b;
                else
                    error('Incorrect size for weight vector assignment');
                end
            case {'X','U'}
                if all(size(b) == size(t.u))
                   t.u = b;
                else
                    error('Incorrect size for factor matrix assignment');
                end
            otherwise
                error(['Field not writable or does not exist: ', s(1).subs]);
        end
    case '()'
        error('Cannot change individual entries in a symktensor.')
    otherwise
        error('Invalid subsasgn.');
end


