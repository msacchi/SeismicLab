function update_tbxref_dir(dirname,varargin)
%UPDATE_TBXREF_DIR Update toolbox link at bottom of each help blurb.
% 
 
%% Parse inputs
params = inputParser;
params.addParameter('Verbose',false);
params.addParameter('Recurse',true);
params.parse(varargin{:});
verbose = params.Results.Verbose;
recurse = params.Results.Recurse;
%% Find all the toolbox files
D = dir(dirname);
if (numel(D) == 0)
    error('ERROR: Cannot find directory!');
end


%% Find all the relevant m-file
if verbose
    fprintf('\n------ Evaluating Directory %s ------\n\n', dirname);
end
for i = 1:numel(D)
    fname = D(i).name;
    ffname = fullfile(dirname,fname);
    if ismember(fname, {'.','..','.git','.gitignore','.gitlab-ci.yml','doc','maintenance','tests'})
        if verbose
            fprintf('Ignoring: %s\n', ffname);
        end
    elseif D(i).isdir
        fprintf('Directory: %s\n', ffname)
        if recurse
            update_tbxref_dir(fullfile(dirname,fname));
        end
    elseif isempty(regexp(fname,'.*\.m$','once'))
        if verbose
            fprintf('Not an m-file: %s\n', ffname)
        end
%     elseif regexp(fname,'tt_','once')
%         if verbose
%             fprintf('Skipping TT file: %s\n', ffname);
%         end
%         update_tbxref_file(ffname)
    else
        if verbose
            fprintf('M-file: %s\n', ffname)
        end
        update_tbxref_file(ffname)
        
    end
end

