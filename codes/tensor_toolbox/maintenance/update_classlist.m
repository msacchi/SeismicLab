function update_classlist(varargin)
%UPDATE_CLASSLIST Updates constructor m-file with list of classes.
%
%   This method assumes there are two lines at the top of the file that
%   define the class. It keeps those in tact. It then searchers for the
%   first line that hase the string "a href" in it, and it replaces all the
%   lines inbetween with the list of functions for the class.
%
%   UPDATE_CLASSLIST('Class',CLASSNAME) updates the filelist at the top of
%   the class constructor file, which is used as documentation for the
%   class when typing 'HELP CLASSNAME'.
%
%   UPDATE_CLASSLIST('Class',CLASSNAME,'Debug',true) doesn't actually
%   overwrite the constructor file. Instead, it just creates the file
%   tmp_CLASSNAME.m in the class directory.
%
%   UPDATE_CLASSLIST updates every class.
%
%   See also CREATE_DIRCONTENTS, CREATE_TOPCONTENTS.
%
%MATLAB Tensor Toolbox. Copyright 2017, Sandia Corporation.


%% Parse inputs
params = inputParser;
params.addParameter('Class',[]);
params.addParameter('Debug',false);
params.addParameter('Copyright',false);
params.parse(varargin{:});

classlist = params.Results.Class;
if isempty(classlist)
    classlist = {'tensor','sptensor','ttensor','ktensor','tenmat',...
        'sptenmat','sumtensor','symtensor','symktensor'};
elseif ~iscell(classlist)
    classlist = {classlist};
end
   
debug = params.Results.Debug;

ttbdir = getfield(what('tensor_toolbox'),'path');

%%
for j = 1:numel(classlist)
    
    % Extract contents
    classname = classlist{j}; 
    classdir = fullfile(ttbdir,strcat('@',classname));
    fname = fullfile(classdir,strcat(classname,'.m'));
    fnametmp = fullfile(classdir,strcat('tmp_',classname,'.m'));
    
    % Extract directory contents
    C = create_dircontents(classdir,'Copyright',params.Results.Copyright);
    
    % Write to main class file
    fprintf('Replacing list of functions in file %s\n', fname);
    fidold = fopen(fname, 'r');
    fidnew = fopen(fnametmp, 'w');

    % Copy first two lines, which are just the class name and
    % description, plus a blank line.
    for i = 1:2
        oldline = fgetl(fidold);
        fprintf(fidnew, '%s\n', oldline);
    end
        
    % Insert new contents into temporary file.
    fprintf(fidnew,'%%%s Methods:\n',upper(classname));
    for i = 1:numel(C)
        fprintf(fidnew,'%%   %s\n',C{i});
    end
    fprintf(fidnew,'%%\n');
        
    % Skip until href found
    while 1
        oldline = fgetl(fidold);
        if ~ischar(oldline)
            error('Never found href line')
        end
        if contains(oldline,'a href=')
            break;
        end
    end
            
    % Just copy everything else
    while 1
        fprintf(fidnew, '%s\n', oldline);
        oldline = fgetl(fidold);
        if ~ischar(oldline), break, end
    end
    
    fclose(fidold);
    fclose(fidnew);
    
    if ~debug
        [s,m] = movefile(fnametmp,fname);
        if (s == 0)
            fprintf('Error renaming file: %s\n',m);
        end
    end
end
 
%


