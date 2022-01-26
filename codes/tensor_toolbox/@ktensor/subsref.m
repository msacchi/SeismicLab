function a = subsref(t,s)
%SUBSREF Subscripted reference for a ktensor.
%
%   Subscripted reference is used to query the components of a ktensor.
%
%   Examples
%   X = ktensor([3; 2], rand(4,2), rand(5,2), rand(3,2));
%   X.lambda %<--returns the lambda array ([3;2]).
%   X.U %<--returns a cell array of 3 matrices.
%   X.U{1} %<--returns the matrix corresponding to the first mode.
%   X(2,3,1) %<--calculates and returns that single element of X.
%
%   See also KTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



switch s(1).type    
    case '.'
        switch s(1).subs
            case 'lambda'
                a = tt_subsubsref(t.lambda,s);
            case {'u','U','a','A'}
                a = tt_subsubsref(t.u,s);
            otherwise
                error(['No such field: ', s(1).subs]);
        end
    case '()'
        if length(s.subs) == 1
            subs = s.subs{1};
            r = length(t.lambda);
            ns = size(subs,1);
            d = size(subs,2);
            b = repmat(reshape(t.lambda,1,r),ns,1);
            for k = 1:d
                tmp = b;
                b = tmp .* t.u{k}(subs(:,k),:);
            end
            a = sum(b,2);
        else
            a = 0;
            for k = 1 : length(t.lambda)
                b = t.lambda(k);
                for i = 1 : length(s.subs)
                    b = b * t.u{i}(s.subs{i},k);
                end
                a  = a + b;
            end
        end
    case '{}'
        a = subsref(t.u,s);
    otherwise
        error('Invalid subsref');
end
