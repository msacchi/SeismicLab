function z = lt(x,y)
%LT Less than for sptensors.
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%% Observations for sparse matrix case.
% The result of a < 5 is sparse.
% The result of a < 0 is sparse.
% The result of a < full(a) is sparse.

%% Case 1: One argument is a scalar
if isscalar(y)
    subs1 = x.subs((x.vals < y),:);
    if y > 0
        subs2 = setdiff(allsubs(x),x.subs,'rows');
    else
        subs2 = [];
    end
    z = sptensor([subs1;subs2],true,size(x));
    return;
end

% Call back with the arguments reversed.
if isscalar(x)
    z = gt(y,x);
    return;
end

%% Case 2: Both x and y are tensors or some sort
% Check that the sizes match
if ~isequal(x.size,y.size)
    error('Size mismatch');
end

% Case 2a: Two sparse tensors
if isa(x,'sptensor') && isa(y,'sptensor')

    % x not zero, y zero
    subs1 = setdiff(x.subs,y.subs,'rows');
    subs1 = subs1(extract(x,subs1) < 0, :);

    % x zero, y not zero
    subs2 = setdiff(y.subs,x.subs,'rows');
    subs2 = subs2(extract(y,subs2) > 0, :);

    % x and y not zero
    subs3 = intersect(x.subs,y.subs,'rows');
    subs3 = subs3(extract(x,subs3) < extract(y,subs3),:);
    
    % assemble
    z = sptensor([subs1;subs2;subs3],true,size(x));
    return;

end

% Case 2b: y is a dense tensor
if isa(y,'tensor')

    % x zero and y > 0
    subs1 = find(y > 0);
    if ~isemtpy(subs1)
        subs1 = setdiff(subs1,x.subs,'rows');
    end
    
    % x and y nonzero
    subs2 = x.subs(x.vals < y(x.subs,'extract'),:);
    
    % assemble
    z = sptensor([subs1;subs2],true,size(x));
    
    return;
    
end

%% Otherwise
error('The arguments must be two sptensors or an sptensor and a scalar.');
