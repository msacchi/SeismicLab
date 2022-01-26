function update_copyright(fname,varargin)
%UPDATE_COPYRIGHT Update pre-3.0 copyright to new version.

%% Setup

%% Parse inputs
params = inputParser;
params.addParameter('Verbose',true);
params.addParameter('Debug',true);
params.parse(varargin{:});
verbose = params.Results.Verbose;
debug = params.Results.Debug;

%% Set up files
[pathstr,name,ext] = fileparts(fname);
fnametmp = fullfile(pathstr, [name '_tmp' ext]);
if (verbose)
    fprintf('Replacing copyright in file %s\n', fname);
end

%% Open files
fidold = fopen(fname, 'r');
fidnew = fopen(fnametmp, 'w');

%% Copy over and replace copyright
while 1
    tline = fgetl(fidold);
    if ~ischar(tline), break, end

    if regexp(tline,'.MATLAB Tensor Toolbox.\w*$')
        fprintf(fidnew,'%%MATLAB Tensor Toolbox. Copyright 2017, Sandia Corporation.\n');
        fprintf(fidnew,'%%https://gitlab.com/tensors/tensor_toolbox.\n');
        continue;
    end
    if regexp(tline,'Copyright 2015, Sandia Corporation'), continue, end
    if regexp(tline,'This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.\w*')
        for i = 1:6
            tline = fgetl(fidold);
        end
        continue;
    end
    fprintf(fidnew, '%s\n', tline);

end

%% Close files
fclose(fidold);
fclose(fidnew);

