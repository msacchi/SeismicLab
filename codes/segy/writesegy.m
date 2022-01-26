function   writesegy(filename,D,H);
%WRITESEGY: Write a SU-SEGY file.
% 
%  Writes Header and Data to a file
%          
%  writesegy(filename,D,H);
%
%  IN  filename: name of file to create
%      D:        data in a matrix  D(time,trace)
%      H:        Header structure containing SEGY headers
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.

 
   FID = fopen(filename,'w','b');      % l,b little or big endian 

   [ns,ntraces] = size(D);             % dimension of D
                                       % ns number of samples per trace

 for k=1:ntraces,
   
   fwrite(FID,H(k).tracl,'int');    
   fwrite(FID,H(k).tracr,'int');     
   fwrite(FID,H(k).fldr,'int');     
   fwrite(FID,H(k).tracf,'int'); 
   fwrite(FID,H(k).ep,'int');       
   fwrite(FID,H(k).cdp,'int');     
   fwrite(FID,H(k).cdpt,'int');    
   fwrite(FID,H(k).trid,'short'); 
   fwrite(FID,H(k).nva,'short');    
   fwrite(FID,H(k).nhs,'short');   
   fwrite(FID,H(k).duse,'short');   
   fwrite(FID,H(k).offset,'int');   
   fwrite(FID,H(k).gelev,'int');  
   fwrite(FID,H(k).selev,'int');    
   fwrite(FID,H(k).sdepth,'int');   
   fwrite(FID,H(k).gdel,'int');     
   fwrite(FID,H(k).sdel,'int');     
   fwrite(FID,H(k).swdep,'int');   
   fwrite(FID,H(k).gwdep,'int'); 
   fwrite(FID,H(k).scalel,'short'); 
   fwrite(FID,H(k).scalco,'short'); 
   fwrite(FID,H(k).sx,'int');       
   fwrite(FID,H(k).sy,'int');       
   fwrite(FID,H(k).gx,'int');       
   fwrite(FID,H(k).gy,'int');       
   fwrite(FID,H(k).counit,'short'); 
   fwrite(FID,H(k).wevel,'short');  
   fwrite(FID,H(k).swevel,'short'); 
   fwrite(FID,H(k).sut,'short');    
   fwrite(FID,H(k).gut,'short');    
   fwrite(FID,H(k).sstat,'short');  
   fwrite(FID,H(k).gstat,'short');  
   fwrite(FID,H(k).tstat,'short');  
   fwrite(FID,H(k).laga,'short');   
   fwrite(FID,H(k).lagb,'short');   
   fwrite(FID,H(k).delrt,'short');  
   fwrite(FID,H(k).muts,'short');   
   fwrite(FID,H(k).mute,'short');   
   fwrite(FID,H(k).ns,'short');  
   fwrite(FID,H(k).dt,'short');  
   fwrite(FID,H(k).gain,'short');   
   fwrite(FID,H(k).igc,'short');    
   fwrite(FID,H(k).igi,'short');   
   fwrite(FID,H(k).corr,'short');   
   fwrite(FID,H(k).sfs,'short');    
   fwrite(FID,H(k).sfe,'short');    
   fwrite(FID,H(k).slen,'short');  
   fwrite(FID,H(k).styp,'short');  
   fwrite(FID,H(k).stas,'short');   
   fwrite(FID,H(k).stae,'short');  
   fwrite(FID,H(k).tatyp,'short');  
   fwrite(FID,H(k).afilf,'short');  
   fwrite(FID,H(k).afils,'short');  
   fwrite(FID,H(k).nofilf,'short'); 
   fwrite(FID,H(k).nofils,'short'); 
   fwrite(FID,H(k).lcf,'short');    
   fwrite(FID,H(k).hcf,'short');   
   fwrite(FID,H(k).lcs,'short');   
   fwrite(FID,H(k).hcs,'short');   
   fwrite(FID,H(k).year,'short');   
   fwrite(FID,H(k).day,'short');    
   fwrite(FID,H(k).hour,'short');   
   fwrite(FID,H(k).minute,'short'); 
   fwrite(FID,H(k).sec,'short');    
   fwrite(FID,H(k).timbas,'short'); 
   fwrite(FID,H(k).trwf,'short');   
   fwrite(FID,H(k).grnors,'short'); 
   fwrite(FID,H(k).grnofr,'short'); 
   fwrite(FID,H(k).grnlof,'short'); 
   fwrite(FID,H(k).gaps,'short');   
   fwrite(FID,H(k).otrav,'short');  
   fwrite(FID,H(k).d1,'float');     
   fwrite(FID,H(k).f1,'float');     
   fwrite(FID,H(k).d2,'float');     
   fwrite(FID,H(k).f2,'float');     
   fwrite(FID,H(k).ungpow,'float'); 
   fwrite(FID,H(k).unscale,'float');
   fwrite(FID,H(k).ntr,'int');     
   fwrite(FID,H(k).mark,'short');  
   fwrite(FID,H(k).unass,'short');
   fwrite(FID,D(:,k) ,'float');
 end

[message,errnum] = ferror(FID);
fclose(FID);
