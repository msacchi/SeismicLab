%% Script to print capability chart for all the tensor toolbox classes.
% Note that we exclude the tensor-as-matrix style classes.

%% Manual class names
    
nclasses = 7;
classnames = {'tensor'; 'sptensor'; 'symtensor'; 'ttensor'; 'ktensor'; 'symktensor'; 'sumtensor'};


%% Get directory contents for each class (omitting constructor)
classmembers = cell(nclasses,1);
functionnames = {};
for i = 1:nclasses
    C = dir(['../@' classnames{i}]);
    nc = length(C);
    tf_tmp = false(nc,1);
    for j = 1:nc
        if (length(C(j).name) > 3) && strcmp(C(j).name(end-1:end),'.m') ...
                && ~strcmp(C(j).name(1:end-2),classnames{i})
            tf_tmp(j) = true;
        end
    end
    C = C(tf_tmp);
    C = arrayfun(@(x) x.name(1:end-2), C, 'UniformOutput', false);
    classmembers{i} = C;
    functionnames = union(functionnames,C);
end

% get membership array
nfunctions = length(functionnames);
tf = false(nfunctions, nclasses);
for i = 1:nclasses
    tf(:,i) = ismember(functionnames, classmembers{i});
end

%% Print out results
cnl = cellfun(@length,classnames);

for i = 1:nfunctions
    if mod(i,20)==1
        fprintf('function     ');
        for i = 1:nclasses
            fprintf('%s ', classnames{i});
        end
        fprintf('\n');
    end
    fprintf('%-12s', functionnames{i});
    for j = 1:nclasses
        for k = 1:4
            fprintf(' ');
        end
        if tf(i,j)
            fprintf('X');
        else
            fprintf('-');
        end
        for k = 5:cnl(j)
            fprintf(' ');
        end
    end
    fprintf('\n');
end
        
        

