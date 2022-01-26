function Y = symmetrize(X,grps,ver)
%SYMMETRIZE Symmetrize a tensor X in specified modes.
%
%   Y = symmetrize(X) will symmetrize a tensor X with respect to all
%   modes so that Y is symmetric with respect to any permutation of
%   indices. The dimensions of Y must be equal accross all modes. The
%   resulting symmetrized tensor is formed by computing the average over
%   all elements in a permutation class.
%
%   Y = symmetrize(X,MODES) will symmetrize a tensor X with respect to the
%   modes specified by the vector MODES of mode indices. The second
%   argument may alternatively be a cell array of vectors of modes for
%   symmetrization.
%
%   NOTE: It is *the same or less* work to just call X = symmetrize(X) then
%   to first check if X is symmetric and then symmetrize it, even if X is
%   already symmetric.
%
%   Examples
%   W = tensor(rand(2,2,2));
%   X = tensor(rand(4,3,3,4));
%   symmetrize(W)
%   symmetrize(X,{[1 4], [2,3]}) %<--Symmetrize in modes [1 4] then [2 3]
%
%   See also TENSOR, TENSOR/ISSYMMETRIC.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


n = ndims(X);
sz = size(X);

%ver is an optional argument specifying the version to use.
if ~exist('ver', 'var')
    ver=0; %By default use new, faster version of issymmetric
end

% Check that grps exists; if not, create it.
if ~exist('grps','var')
    grps = 1:n;
end

% Check that grps is a cell array.
if ~iscell(grps)
    grps = {grps};
end

switch ver

    case 0 % New version (default!)
        
        ngrps = length(grps);
        for i = 1:ngrps

            thisgrp = grps{i};
            
            % Check tensor dimensions for compatibility with symmetrization
            if ~all( sz(thisgrp(1)) == sz(thisgrp) )
                error('TTB:Tensor:BadModes','Dimension mismatch for symmetrization');
            end
            
            % Check for no overlap in the sets
            if i < ngrps  
                if ~all(isempty(intersect(thisgrp,grps{i+1:end})))
                    error('TTB:Tensor:BadModes','Cannot have overlapping symmetries');
                end
            end            
            
            % Construct matrix ind where each row is the multi-index for
            % one element of X
            idx = tt_ind2sub(size(X), (1:numel(X.data))');
            
            % Find reference index for every element in the tensor - this
            % is to its index in the symmetrized tensor. This puts every
            % element into a 'class' of entries that will be the same under
            % symmetry.
            classidx = idx;
            classidx(:,thisgrp) = sort(idx(:,thisgrp),2);
            linclassidx = tt_sub2ind(size(X), classidx);
            
            % Skip if its already symmetric
            if all(X.data(:) == X.data(linclassidx));
                continue;                
            end
             
            % Take average over all elements in the same class
            classsum = accumarray(linclassidx, X.data(:));
            classnum = accumarray(linclassidx, 1);            
            avg = classsum ./ classnum;              
        
            % Fill in each entry with its new symmetric version 
            newdata = avg(linclassidx);  
            X.data = reshape(newdata,[size(X) 1]);
        end
        
        % Final result
        Y = X;
        
    case 1  % The original version of the algorithm
        
        % Check tensor dimensions for compatibility with symmetrization
        ngrps = length(grps);
        for i = 1:ngrps
            dims = grps{i};
            for j = dims(2:end)
                if sz(j) ~= sz(dims(1))
                    error('Dimension mismatch for symmetrization');
                end
            end
        end
        
        % Check for no overlap in the sets
        for i = 1:ngrps
            for j = i+1:ngrps
                if ~isempty(intersect(grps{i},grps{j}))
                    error('Cannot haver overlapping symmetries');
                end
            end
        end
        
        % Create the combinations for each symmetrized subset
        combos = cell(ngrps,1);
        for i = 1:ngrps
            combos{i} = perms(grps{i});
        end
        
        % Create all the permuations to be averaged
        total_perms = prod(cellfun(@length,combos));
        sym_perms = repmat(1:n, total_perms, 1);
        for i = 1:ngrps
            ntimes = prod(cellfun(@length,combos(1:i-1)));
            ncopies = prod(cellfun(@length,combos(i+1:end)));
            nelems = length(combos{i});
            
            idx = 1;
            for j = 1:ntimes
                for k = 1:nelems
                    for l = 1:ncopies
                        sym_perms(idx,grps{i}) = combos{i}(k,:);
                        idx = idx + 1;
                    end
                end
            end
        end
        
        % Create an average tensor
        Y = tenzeros(size(X));
        for i = 1:total_perms
            Y = Y + permute(X,sym_perms(i,:));
        end
        Y = Y / total_perms;
        
        % It's not *exactly* symmetric due to oddities in differently ordered
        % summations and so on, so let's fix that.
        % Idea borrowed from Gergana Bounova:
        % http://www.mit.edu/~gerganaa/downloads/matlab/symmetrize.m
        for i = 1:total_perms
            Z = permute(Y,sym_perms(i,:));
            Y.data(:) = max(Y.data(:),Z.data(:));
        end
        
    otherwise
        error('incorrect version specification');
end



