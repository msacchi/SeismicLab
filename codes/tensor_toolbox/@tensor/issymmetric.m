function [tf,all_diffs,all_perms] = issymmetric(X,grps,ver)
%ISSYMMETRIC Verify that a tensor X is symmetric in specified modes.
%
%   TF = ISSYMMETRIC(X) returns true if X is exactly symmetric for every
%   permutation of its modes.
%
%   [TF,DIFFS,PERMS] = ISSYMMETRIC(X) also returns that maximum difference
%   in DIFFS for each permutation in PERMS (one permutation per row).
%
%   [...] = ISSYMMETRIC(X,IDX) checks symmetry with respect to the modes
%   specified in IDX, which can be an array of indices or a cell array of
%   arrays of symmetric indices.
%
%   Examples
%   W = tensor(rand(3,3,3));
%   issymmetric(W,[1 2]) %<--Checks for symmetry in modes [1 2]. False here
%
%   See also TENSOR, SYMMETRIZE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


n = ndims(X);
sz = size(X);

%ver is an optional argument specifying the version to use.
if ~exist('ver', 'var')
    ver = 0; %By default use new, faster version of issymmetric
end

% Check that grps exists; if not, create it.
if ~exist('grps','var')
    grps = 1:n;
end

if nargout > 1
    ver = 1; % User requested permutation and difference information.
end

% Check that grps is a cell array.
if ~iscell(grps)
    grps = {grps};
end

%Substantially different routines are called depending on whether the user
%requests the permutation information. If permutation is required (or requested)
%the algorithm is much slower
switch ver
    
    case 0 %use new algorithm
    
        for i = 1:length(grps)
            
            % Extract current group
            thisgrp = grps{i};
            
            % Check tensor dimensions first
            if ~all( sz(thisgrp(1)) == sz(thisgrp) )
                tf = false;
                return;
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
            
            % Compare each element to its class exemplar
            if any(X.data(:) ~= X.data(linclassidx));
                tf = false;
                return
            end
        end
        
        % We made it past all the tests!
        tf = true;
        return
        
    case 1 %Use older algorithm
        
        % Check tensor dimensions for compatibility with symmetrization
        for i = 1:length(grps)
            dims = grps{i};
            for j = dims(2:end)
                if sz(j) ~= sz(dims(1))
                    tf = false;
                    return;
                end
            end
        end
        
        % Check actual symmetry.
        cnt = sum(cellfun(@(x) factorial(length(x)), grps));
        all_diffs = zeros(cnt,1);
        all_perms = zeros(cnt,n);
        idx = 1;
        for i = 1:length(grps)
            
            % Compute the permutations for this group of symmetries
            p = perms(grps{i});
            
            for j = 1:size(p,1)
                
                % Create the permutation to check
                q = 1:n;
                q(grps{i}) = p(j,:);
                
                % Save the permutation
                all_perms(idx,:) = q;
                
                % Do the permutation and see if it's a match. If it's not a match,
                % record the difference.
                Y = permute(X,q);
                if isequal(X.data,Y.data)
                    all_diffs(idx) = 0;
                else
                    all_diffs(idx) = max(abs(X.data(:)-Y.data(:)));
                end
                
                % Increment the index
                idx = idx + 1;
                
            end
            
        end
        
    otherwise
        error('Incorrect version specification');
end
tf = all(all_diffs == 0);
end
