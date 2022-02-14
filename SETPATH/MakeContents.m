function MakeContents(aString, flags)
% MAKECONTENTS makes Contents file in current working directory.
%   MAKECONTENTS(STRING,FLAGS) creates a standard "Contents.m" file in the
%   current directory by assembling the first comment (H1) line in
%   each function found in the current working directory.  If a 
%   "Contents.m" file exists, it is renamed to "Contents.old", before
%   a new "Contents.m" file is created.  STRING is inserted as the 
%   first line of the "Contents.m" file;  if omitted, a blank line 
%   is inserted instead.  The function changes permission on the 
%   resultant "Contents.m" file to rw-r--r-- on Unix systems.
%
%   FLAGS can contain none or more of
%      'n'    - suppress path name in Contents file
%      'f'    - include first word of first line (excluded by default)
%      'F'    - include first word only if not capitalized
%      'c'    - use filename 'contents.m' instead of 'Contents.m'
%      'r'    - recursively list subdirectory contents also
%
% Updated 29 June 2000.
% Revised to recurse down directories, handle options by
% Matthew Brett; 28 June 2003
%
% See also CONTENTS.

% Author(s): L. Bertuccioli 
%            A. Prasad

% Based on mkcontents.m by Denis Gilbert

% Default value of input string
if nargin < 1,
     aString =' ';
end
if nargin < 2
  flags = '';
end
if isempty(flags)
  flags = ' ';
end

cont_file = 'Contents.m';
if any(flags == 'c'), cont_file = lower(cont_file); end

disp(['Creating "' cont_file '" in ' pwd])
if exist([pwd filesep cont_file]) ~= 0 
     copyfile(cont_file,[cont_file(1:end-1) 'old']);
     delete(cont_file)
end

% Header lines
line1 = ['% ' aString];
fcontents = fopen(cont_file,'w'); 
fprintf(fcontents,'%s\n',line1);     
if ~any(flags == 'n') % add path to directory
  line2 = ['% Path ->  ' pwd];
  fprintf(fcontents,'%s\n',line2);     
end

% do write
do_list('.', fcontents, flags);
fclose(fcontents);

% If Unix. change permissions on Contents.m file
if isunix
  unix(['chmod go+r ' cont_file]);
end
return

function do_list(dirname, fcontents, flags);
% subfunction to print contents information into contents file

% Parse directory for mfile names and directories
dirlist = dir(dirname);
dir_i = [dirlist.isdir];
dirnames = {dirlist(dir_i).name};
% remove . and .. from directory list
dirnames = dirnames(~(strcmp('.', dirnames) | strcmp('..', dirnames)));
filenames = {dirlist(~dir_i).name};
mfiles = {};
% Find .m files, excluding Contents file
for f = 1:length(filenames)
  fname = filenames{f};
  if length(fname) > 2 & ...
	strcmp(fname(end-1:end), '.m') & ...
	       ~strcmpi('contents.m', fname)
    mfiles = [mfiles {fname}];
  end
end

fprintf(fcontents,'%%\n'); 

if strcmp(dirname, '.')
  dirlab = dirname(2:end);
else
  dirlab = [dirname(3:end) filesep];
end
maxlen = size(char(mfiles),2) + length(dirlab);

% Write first lines to Contents.m if they exist
for i = 1:length(mfiles)
   mfile = mfiles{i};
   fid=fopen(fullfile(dirname, mfile),'r');
   if fid == -1
      error(['Could not open file: ' fullfile(dirname, mfile)]);
   end
   aLine = '';
   while(isempty(aLine))
     aLine = fgetl(fid);
   end
   if length(aLine) > 7  & strcmp(aLine(1:8),'function') == 1,
	count_percent = 0;
	while count_percent < 1 & feof(fid)==0; 
	     line = fgetl(fid);
	     if length(line) > 0 
		  if ~isempty(findstr(line,'%')) 
		       count_percent = count_percent + 1;
		       rr=line(2:length(line));
		       % get first word
		       [tt,rr2]=strtok(line(2:length(line)));
		       % ? remove first word
		       if ~(any(flags == 'f') | ...
			 (any(flags == 'F') & strcmp(lower(tt),tt)))
			 rr = rr2;
		       end
		       rr = fliplr(deblank(fliplr(rr)));
		       fn = [dirlab strtok(mfile,'.')];
		       n = maxlen - length(fn) - 1;
		       line = ['%   ' fn blanks(n) '- ' rr];
		       fprintf(fcontents,'%s\n',line);
		  end % if ~isempty
	     end % if length
	     if feof(fid)==1  
		  fn = [dirlab strtok(mfile,'.')];
		  n = maxlen - length(fn) - 1;
		  line = ['%   ' fn blanks(n) '- (No help available)'];
		  fprintf(fcontents,'%s\n',line); 
	     end % if feof
	end % while
   end % if strcmp
   fclose(fid);
end
% recurse down directory tree if required
if any(flags == 'r')
  for d = 1:length(dirnames)
    do_list(fullfile(dirname, dirnames{d}), fcontents, flags)
  end
end
return
