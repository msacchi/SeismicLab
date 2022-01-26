function t = subsasgn(t,s,b)
%SUBSASGN Subscripted assignment for ktensor.
%
%   Subscripted assignment can be used to alter the lambda vector or the
%   factor matrices of a ktensor. The entire factor matrix or weight vector
%   must be provided.
%
%   Examples
%   X = ktensor(ones(4,1), rand(2,4), rand(3,4), rand(4,4));
%   X.lambda = 2*ones(4,1) %<--Redefine weight vector
%   X.U{1} = zeros(2,4) %<--Redefine first factor matrix
%   X.U = {zeros(2,4), zeros(3,4), zeros(4,4)} %<--Redefine factor matrices
%
%   See also KTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



switch s(1).type
    case '.'
        switch s(1).subs
            case 'lambda'
                if length(s) == 1
                    t = ktensor(b, t.u);
                else
                    newlambda = subsasgn(t.lambda, s(2:end), b);
                    t = ktensor(newlambda, t.u);
                end
            case {'u','U'}
                if length(s) == 1
                    t = ktensor(t.lambda, b);
                else
                    tmpu = subsasgn(t.u, s(2:end), b);
                    t = ktensor(t.lambda, tmpu);
                end
            otherwise
                error(['No such field: ', s(1).subs]);
        end
    case '()'
        error('Cannot change individual entries in a ktensor.')
    case '{}'
        new_s(1).type = '.';
        new_s(1).subs = 'u';
        new_s(2:length(s)+1) = s;
        t = subsasgn(t, new_s, b);
    otherwise
        error('Invalid subsasgn.');
end


