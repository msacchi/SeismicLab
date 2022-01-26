function B = subsref(A, S)
%SUBSREF Subreference function for symtensor.
%
%   Examples:
%   X = symtensor(1:20,3,4);
%   X.val %<- Returns the distinct values
%   X.m %<- Tensor order
%   X.n %<- Tensor dimension
%   X(5) %<- Fifth distinct element
%   X((1:4)') %<- Linear indexing of distinct values
%   X(1,2,1) %<- Returns X(1,2,1) = X(1,1,2) = 2nd element
%   X([1 1 2;3 2 1]) %<- Return two elements
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

sz=size(S);

switch S(1).type
    case '{}'
        error('Cell contents reference from a non-cell array object.')
    case '.'
        fieldname = S(1).subs;
        switch fieldname
            case {'val','vals'}
                if sz(2) == 2  % Can query into the val array
                    B = subsref(A.val, S(2));
                else
                    B = A.val;
                end
                
            case 'm'
                B = A.m;  % Query the number of modes
            case 'n'
                B = A.n;  % Query the dimension
            otherwise
                error(['No such field in symtensor: ', fieldname]);
        end
    case '()'
        if (A.m > 1) && (numel(S.subs) == A.m)
            if ~all(cellfun(@(x) isscalar(x) && isnumeric(x), S.subs))
                error('Invalid indexing for symktensor');
            end
            newS = S;
            newS.subs = cell(1,1);
            newS.subs{1} = cell2mat(S.subs);
            B = subsref(A,newS);
            return;
        elseif  numel(S.subs) == 1 
            if size(S.subs{1},2) == 1 % Linear indexing (into value array)
                B = subsref(A.val,S);
                return
            end
            if size(S.subs{1},2) == A.m % Query is a matrix whose rows are indices
                qsubs = sort(S.subs{1},2); % Sort the indices
                asubs = indices(A);
                [tf,loca] = ismember(qsubs,asubs,'rows');
                if ~all(tf)
                    error('TTB:Symtensor:BadIdx', 'Invalid Indexing');
                end
                B = A.val(loca);
                return
            end
            error('Invalid Indexing')
        else
            error('Invalid indexing');
        end
    otherwise
        error('Invalid indexing for symktensor');
end


end

