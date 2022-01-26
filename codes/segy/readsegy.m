function [D,H] = readsegy(filename,hw,min,max)
%READSEGY: Read SU-SEGY data.
%          The data and headers can be  extracted  in 
%          a range given by hw (header words).
%   
%  [D,H] = readsegy(filename,hw,min,max) returns the data D and
%
%  IN   filename: name of segy
%       hw:       header word to limit number of traces to read
%       min, max: min and maximum value of hw
%
%  OUT  D:        the data (in a matrix) 
%       H:        the header in a structure
%
%
%  Examples:    
%
%    [D,H] = readsegy('data_with_noise.su');
%    [nt,nh] = size(D);
%    offset = [H.offset];
%    dtsec = H(1).dt/1000/1000;
%    imagesc(offset,[0:1:nt-1]*dtsec,D);
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.

 
  FID = fopen(filename,'r','b');       % l: little, b: big endian

   segy=segy_struct;                   % load the definitions of 
                                       % the header words

   count=count_struct;                 % load the position of 
                                       % each word in the header (in bytes)

   status = fseek(FID,count.ns,'bof'); % go to the beggining of file

   ns = fread(FID,1,segy.ns);          % read ns from first trace
                                       % ns is the number of traces per trace

   total = 60+ns;                      % total nuber of 4bytes words


   max_traces=9999999;                 % maximum number of traces (will
                                       % stop before). The variable status
                                       % will make the code stop when
                                       % the end of file is reached  

 if nargin~=1;
   hc=eval(strcat('count.',hw));       % assigned the header word required  
   hp=eval(strcat('segy.',hw));        % to extract the traces.
 j = 1;                                % counter 
 for k =1:max_traces
   position = total*(k-1)*4+hc;        % where in bytes is the header word
   status = fseek(FID,position,'bof'); 
  if status == 0                       % stop when end of file is reached
    w = fread(FID,1,hp);
     if  w>=min                        % pick traces with a given range
      if w<=max                        % of the desired header word
       position = total*(k-1)*4+count.trace;
       status = fseek(FID,position,'bof'); 
       trace = fread(FID,ns,segy.trace);
       j = j + 1;  
       D(:,j-1)  = trace(:);           % load traces into D
       H(j-1)  = header(FID,ns,k);     % load each header in a structure H
  if  rem(j-1,1000)==0;  
       fprintf('Traces = % d',j-1)
  end
      end
     end
   else
 fprintf(' Number of traces =  %d',j-1)
  return
 end
 end 


 else 

% when no hw and bounds are given, reads evrything

 for k =1:max_traces
       position = total*(k-1)*4+count.trace;
       status = fseek(FID,position,'bof'); 
       if status == 0                 
       trace = fread(FID,ns,segy.trace);
       D(:,k)  = trace(:);           % load traces into D
       H(k)  = header(FID,ns,k);     % load each header in a structure H
  if  rem(k,500)==0;  
       fprintf(' Traces = %d\n',k)
  end
   else
   fprintf(' Number of traces =  %d\n',k-1)
  return
  end
 end 
Ntraces=k-1;
 end 



 [message,errnum] = ferror(FID)
 fclose(FID);

