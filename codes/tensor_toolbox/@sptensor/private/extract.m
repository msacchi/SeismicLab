function a = extract(t,srchsubs)
%EXTRACT Extract value for a sptensor. 
%
%   EXTRACT(X,SUBS) returns a list of values.
%
%   See also SPTENSOR/SUBSREF.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>




% Check range of requested subscripts
p = size(srchsubs,1);

% Check that all subscripts are positive and less than the max
invalid = (srchsubs < 0) | (srchsubs > ones(p,1)*t.size);
badloc = find(sum(invalid,2));       
if ~isempty(badloc)
    fprintf('The following subscripts are invalid: \n');
    badsubs = srchsubs(badloc,:);
    badidx = tt_sub2ind(size(t),badsubs);
    for i = 1:numel(badloc)
        fprintf('\tsubscript = %s (linear index = %d)\n',...
                    tt_intvec2str(badsubs(i,:)), badidx(i));
    end
    error('Invalid subscripts');
end

% Set the default answer to zero
a = zeros(p,1);

% Find which indices already exist and their locations
[tf,loc] = ismember(srchsubs,t.subs,'rows');

% Fill in the non-zero elements in the answer
nzsubs = find(tf);
a(nzsubs,1) = t.vals(loc(nzsubs));

