function t = subsasgn(t,s,b)
%SUBSASGN Subscripted assignment for a ttensor.
%   
%   Subscripted assignment can be used to alter the core tensor or the
%   factor matrices in a ttensor. The entire factor matrix or tensor must be
%   provided.
%
%   Examples
%   X = ttensor(tensor(rand(2,2,2)), rand(4,2), rand(5,2), rand(3,2));
%   X.core = tensor(ones(2,2,2)) %<--Change core tensor to all ones
%   X.U{2} = zeros(4,2) %<--Change 2nd factor matrix to zeros
%   X.U = {zeros(4,2),ones(5,2), randn(3,2)} %<--Change all matrices at once
%
%   See also TTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


switch s(1).type
    case '.'
        switch s(1).subs
            case {'core','lambda'}
                if length(s) == 1
                    t = ttensor(b, t.u);
                else
                    tmpcore = subsasgn(t.core, s(2:end), b);
                    t = ttensor(tmpcore, t.u);
                end
            case {'u','U'}
                if length(s) == 1
                    t = ttensor(t.core, b);
                else
                    tmpu = subsasgn(t.u, s(2:end), b); %Refine in U
                    t = ttensor(t.core, tmpu);
                end
            otherwise
                error(['Cannot change field ', s.subs, ' directly.']);
        end
    case '()'
        error('Cannot change individual entries in ttensor.')
    case '{}'
        new_s(1).type = '.';
        new_s(1).subs = 'u';
        new_s(2:length(s)+1) = s;
        t = subsasgn(t, new_s, b);
    otherwise
        error('Invalid subsasgn.');
end


