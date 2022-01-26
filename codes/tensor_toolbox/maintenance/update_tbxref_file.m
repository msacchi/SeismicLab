function update_tbxref_file(fullfname)
%UPDATE_TBXREF_FILE Update the toolbox reference in the given file.




%% Open file
[pathstr,fname,ext] = fileparts(fullfname);
fidold = fopen(fullfname,'r');

%% Old stuff to eliminate
elimstrs = {...
'% This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.',...
'% http://www.sandia.gov/~tgkolda/TensorToolbox.',...
'% Copyright (2015) Sandia Corporation. Under the terms of Contract',...
'% DE-AC04-94AL85000, there is a non-exclusive license for use of this',...
'% work by or on behalf of the U.S. Government. Export of this data may',...
'% require a license from the United States Government.',...
'% The full license terms can be found in the file LICENSE.txt'};

%% New reference
newtbxref = '%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>';


%% Read and process file

% 0 = before function declaration
% 1 = function declaration
% 2 = help title line
% 3 = blank help line
% 4 = indented help contents
% 5 = final help line (not indented) 
% 5.XX = copyright line (extra final help line)
% 6 = after help section
state = 0;
tnum = 0;
fchange = false;
fnum = 0;
fcontents = cell(1000,0);
while 1
    
    % Get next line of text
    tline = fgetl(fidold);
    
    % Check for end of file
    if ~ischar(tline)
        break
    else
        tnum = tnum + 1;
        fnum = fnum + 1;
        fcontents{fnum} = tline;
    end  
   
    if state == 0 
        % This is the state before the function declaration, which occurs
        % in the class defintion files. Once we find the function
        % declaration, we change the state to 1.
        if regexp(tline,'function') == 1
            state = 1;
        end
    elseif (state == 1) 
        % The last line was the function declaration, so the next line
        % *must* be the title help line, which has a specific format. If
        % it's not this line, then we throw an error. Otherwise, we check
        % the formatting of the line and change the state to 2.
        if regexp(tline,'^%.*')
            titleline = regexp(tline,'%([A-Z_0-9]*)\s*(.*)','tokens');
            titlename = lower(titleline{1}{1});
            titledesc = titleline{1}{2};
            
            if ~isequal(titlename,fname)
                mywarning('Filename/description mismatch', fullfname, tnum, tline);
            elseif isempty(regexp(titledesc,'\.\s*$','once'))
                mywarning('Missing final period for description', fullfname, tnum, tline);
            end
        else
            error('function declaration not followed by title line')
        end
        state = 2;
    elseif (state > 1) && (state < 6)
        
        % The prior line was a help line.
        
        % Split up the current line into leading spaces and text
        foo = regexp(tline,'^%(\s*)(.*)','tokens');
        
        if isempty(foo)
            % This is not a help line, signaling the end of the help
            % section. Put up an error if the prior line was not the
            % closing help line (state >= 5). 
            if state < 5
                mywarning('Missing final help section line', fullfname, tnum, tline);
            end
            state = 6;
        else
            % We're in the help section, but after the title help line.
            
            bar = foo{1};
            
            if isempty(bar{1}) && ~isempty(bar{2}) && state < 5
                
                % Encountered the closing help line
                state = 5;
                if ~strcmp(tline,newtbxref)
                    fcontents{fnum} = newtbxref;
                    fchange = true;
                end
           
            elseif isempty(bar{1}) && strncmp(bar{2},'Copyright',9)
                
                % Sometimes there is an extra closing help line that starts
                % with Copyright. This will just be deleted when we write
                % over the closing line.
                state = state + 0.01;
                fnum = fnum - 1 ;
                fchange = true;

            
            else
                
                if state > 5
                    % We've already seen what we thought was the final help
                    % line, so more help lines may indicate a problem.
                    mywarning('Final help line was not final?!?!', fullfname, tnum, tline);
                end
             
               if isempty(bar{2})
                   % Empty help line
                   state = 3;
               else
                   if state == 2 % Prior line was title help line
                       mywarning('Missing blank line following title help line', fullfname, tnum, tline);
                   end
                   if length(bar{1}) >= 3 % Normal help line
                       state = 4;
                   elseif ~isempty(bar{1}) % Not last line
                       mywarning('Too few spaces in help line', fullfname, tnum, tline);
                   else % Last line
                       if state ~= 3 
                           mywarning('Missing blank line before final help line in %s',fullfname, tnum, tline);
                       end
                       state = 5;
                   end
               end
            end
        end
        
    else
        
        if ismember(tline, elimstrs)
            fnum = fnum - 1 ;
            fchange = true;
        end
        
    end
    
end
 
%% Clean up
fclose(fidold);

%% Print contents
if fchange
    fidnew = fopen(fullfname,'w');
    for i = 1:fnum
        fprintf(fidnew,'%s\n', fcontents{i});
    end
    fclose(fidnew);
end

%%
function mywarning(wstr, fname, tnum, tline)

fprintf('*** WARNING: %s in file %s ***\n', wstr, myedit(fname))
if exist('tnum','var')
    fprintf('Line %d: %s\n', tnum, tline);
end

function str = myedit(fname)

str = sprintf('<a href="matlab:edit ''%s''">%s</a>',fname,fname);


