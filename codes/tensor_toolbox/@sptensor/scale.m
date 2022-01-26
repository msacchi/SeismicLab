function t = scale(t,s,dims)
%SCALE Scale along specified dimensions for sparse tensors.
%
%   Y = SCALE(X,S,DIMS) scales the sparse tensor X along the
%   dimension(s) specified in DIMS using the scaling data in S. If
%   DIMS contains only one dimensions, then S can be a column
%   vector. Otherwise, S should be a tensor or sparse tensor.
%
%   Examples
%   X = ones(sptenrand([3 4 5], 10))
%   S = 10 * [1:5]'; Y = scale(X,S,3)
%   S = tensor(10 * [1:5]',5); Y = scale(X,S,3)
%   S = tensor(1:12,[3 4]); Y = scale(X,S,[1 2])
%   S = tensor(1:12,[3 4]); Y = scale(X,S,-3)
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','collapse_scale_doc.html')))">Documentation page for collapsing and scaling tensors</a>
%
%   See also SPTENSOR, SPTENSOR/COLLAPSE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



dims = tt_dimscheck(dims,ndims(t));

switch(class(s))
    case {'tensor'}
        if ~isequal(size(s), t.size(dims))
            error 'Size mismatch';
        end
        t.vals = t.vals .* s(t.subs(:,dims),'extract');
    case {'sptensor'}
        if ~isequal(s.size, t.size(dims))
            error 'Size mismatch';
        end
        t.vals = t.vals .* extract(s,(t.subs(:,dims)));
    case {'double'}
        if size(s,1) ~= t.size(dims)
            error 'Size mismatch';
        end
        t.vals = t.vals .* s(t.subs(:,dims));
    otherwise
        error('Invalid scaling factor');
end
