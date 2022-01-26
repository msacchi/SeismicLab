function [D,H] = ssort(D,H,hw)
%SSORT: Sort data using a header word hw.
%
%  [D,H] = ssort(D,H,hw)
%
%  IN   D,H: Data and Headers in su-segy format
%            (read them  read_segy.m)
%       hw : One of the following words
%                      'cdp'
%                      'offset'
%                      'sx'
%                      'gx'
%                      'tracl'
%                      'tracr'
%
%  OUT  D,H: Data and headers after sorting
%
%  NOTE: This is a poor-man version of what should be an important part
%        of SeismicLab. Don't try to use this code to sort large data
%        sets.
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


%

 if ( nargin < 3 ); hw = 'cdp'; end;

      if strcmp(hw,'cdp');    index =  6; end;
      if strcmp(hw,'offset'); index = 12; end;
      if strcmp(hw,'sx');     index = 22; end;
      if strcmp(hw,'gx');     index = 24; end;
      if strcmp(hw,'sy');     index = 23; end;
      if strcmp(hw,'gy');     index = 25; end;
      if strcmp(hw,'tracl');  index =  1; end;
      if strcmp(hw,'tracr');  index =  2; end;

   [ns,ntraces] = size(D);

   tmp = [H.tracl];                   A(1,:) = tmp(:)';
   tmp = [H.tracr];                   A(2,:) = tmp(:)';  
   tmp = [H.fldr];                    A(3,:) = tmp(:)';    
   tmp = [H.tracf];                   A(4,:) = tmp(:)';
   tmp = [H.ep];                      A(5,:) = tmp(:)';      
   tmp = [H.cdp];                     A(6,:) = tmp(:)';    
   tmp = [H.cdpt];                    A(7,:) = tmp(:)';   
   tmp = [H.trid];                    A(8,:) = tmp(:)';
   tmp = [H.nva];                     A(9,:) = tmp(:)';    
   tmp = [H.nhs];                    A(10,:) = tmp(:)';
   tmp = [H.duse];                   A(11,:) = tmp(:)';
   tmp = [H.offset];                 A(12,:) = tmp(:)';
   tmp = [H.gelev];                  A(13,:) = tmp(:)';
   tmp = [H.selev];                  A(14,:) = tmp(:)';   
   tmp = [H.sdepth];                 A(15,:) = tmp(:)';
   tmp = [H.gdel];                   A(16,:) = tmp(:)';
   tmp = [H.sdel];                   A(17,:) = tmp(:)';    
   tmp = [H.swdep];                  A(18,:) = tmp(:)';  
   tmp = [H.gwdep];                  A(19,:) = tmp(:)';
   tmp = [H.scalel];                 A(20,:) = tmp(:)';
   tmp = [H.scalco];                 A(21,:) = tmp(:)';
   tmp = [H.sx];                     A(22,:) = tmp(:)';
   tmp = [H.sy];                     A(23,:) = tmp(:)';
   tmp = [H.gx];                     A(24,:) = tmp(:)';
   tmp = [H.gy];                     A(25,:) = tmp(:)';      
   tmp = [H.counit];                 A(26,:) = tmp(:)';
   tmp = [H.wevel];                  A(27,:) = tmp(:)';
   tmp = [H.swevel];                 A(28,:) = tmp(:)';
   tmp = [H.sut];                    A(29,:) = tmp(:)';
   tmp = [H.gut];                    A(30,:) = tmp(:)';   
   tmp = [H.sstat];                  A(31,:) = tmp(:)';
   tmp = [H.gstat];                  A(32,:) = tmp(:)';
   tmp = [H.tstat];                  A(33,:) = tmp(:)'; 
   tmp = [H.laga];                   A(34,:) = tmp(:)';
   tmp = [H.lagb];                   A(35,:) = tmp(:)';
   tmp = [H.delrt];                  A(36,:) = tmp(:)'; 
   tmp = [H.muts];                   A(37,:) = tmp(:)';  
   tmp = [H.mute];                   A(38,:) = tmp(:)';
   tmp = [H.ns];                     A(39,:) = tmp(:)';   
   tmp = [H.dt];                     A(40,:) = tmp(:)';  
   tmp = [H.gain];                   A(41,:) = tmp(:)';  
   tmp = [H.igc];                    A(42,:) = tmp(:)';   
   tmp = [H.igi];                    A(43,:) = tmp(:)';  
   tmp = [H.corr];                   A(44,:) = tmp(:)';   
   tmp = [H.sfs];                    A(45,:) = tmp(:)';   
   tmp = [H.sfe];                    A(46,:) = tmp(:)';   
   tmp = [H.slen];                   A(47,:) = tmp(:)'; 
   tmp = [H.styp];                   A(48,:) = tmp(:)';
   tmp = [H.stas];                   A(49,:) = tmp(:)';  
   tmp = [H.stae];                   A(50,:) = tmp(:)';
   tmp = [H.tatyp];                  A(51,:) = tmp(:)';  
   tmp = [H.afilf];                  A(52,:) = tmp(:)'; 
   tmp = [H.afils];                  A(53,:) = tmp(:)'; 
   tmp = [H.nofilf];                 A(54,:) = tmp(:)';
   tmp = [H.nofils];                 A(55,:) = tmp(:)';
   tmp = [H.lcf];                    A(56,:) = tmp(:)';  
   tmp = [H.hcf];                    A(57,:) = tmp(:)';
   tmp = [H.lcs];                    A(58,:) = tmp(:)';  
   tmp = [H.hcs];                    A(59,:) = tmp(:)';
   tmp = [H.year];                   A(60,:) = tmp(:)';  
   tmp = [H.day];                    A(61,:) = tmp(:)';   
   tmp = [H.hour];                   A(62,:) = tmp(:)'; 
   tmp = [H.minute];                 A(63,:) = tmp(:)';
   tmp = [H.sec];                    A(64,:) = tmp(:)';
   tmp = [H.timbas];                 A(65,:) = tmp(:)';
   tmp = [H.trwf];                   A(66,:) = tmp(:)';  
   tmp = [H.grnors];                 A(67,:) = tmp(:)';
   tmp = [H.grnofr];                 A(68,:) = tmp(:)';
   tmp = [H.grnlof];                 A(69,:) = tmp(:)';
   tmp = [H.gaps];                   A(70,:) = tmp(:)';   
   tmp = [H.otrav];                  A(71,:) = tmp(:)';  
   tmp = [H.d1];                     A(72,:) = tmp(:)';     
   tmp = [H.f1];                     A(73,:) = tmp(:)';    
   tmp = [H.d2];                     A(74,:) = tmp(:)';    
   tmp = [H.f2];                     A(75,:) = tmp(:)';    
   tmp = [H.ungpow];                 A(76,:) = tmp(:)';
   tmp = [H.unscale];                A(77,:) = tmp(:)';
   tmp = [H.ntr];                    A(78,:) = tmp(:)';    
   tmp = [H.mark];                   A(79,:) = tmp(:)'; 
%  tmp2 = [H.unass]  not used in this prog.        

% Augmented matrix (data plus headers) 

   AD = [A' D']';

% Default sorting

% do the sorting

  AD = sortrows(AD',index)';

% put back headers into a structure 

 for k=1:ntraces

   H(k).tracl= AD(1,k);
   H(k).tracr= AD(2,k);  
   H(k).fldr=AD(3,k);    
   H(k).tracf=AD(4,k);
   H(k).ep=AD(5,k);
   H(k).cdp=AD(6,k);
   H(k).cdpt=AD(7,k);   
   H(k).trid=AD(8,k);
   H(k).nva=AD(9,k);    
   H(k).nhs=AD(10,k);
   H(k).duse=AD(11,k);
   H(k).offset=AD(12,k);
   H(k).gelev=AD(13,k);
   H(k).selev=AD(14,k);   
   H(k).sdepth=AD(15,k);
   H(k).gdel=AD(16,k);
   H(k).sdel=AD(17,k);    
   H(k).swdep=AD(18,k);  
   H(k).gwdep=AD(19,k);
   H(k).scalel=AD(20,k);
   H(k).scalco=AD(21,k);
   H(k).sx=AD(22,k);
   H(k).sy=AD(23,k);
   H(k).gx=AD(24,k);
   H(k).gy=AD(25,k);      
   H(k).counit=AD(26,k);
   H(k).wevel=AD(27,k);
   H(k).swevel=AD(28,k);
   H(k).sut=AD(29,k);
   H(k).gut=AD(30,k);   
   H(k).sstat=AD(31,k);
   H(k).gstat=AD(32,k);
   H(k).tstat=AD(33,k); 
   H(k).laga=AD(34,k);
   H(k).lagb=AD(35,k);
   H(k).delrt=AD(36,k); 
   H(k).muts=AD(37,k);  
   H(k).mute=AD(38,k);
   H(k).ns=AD(39,k);   
   H(k).dt=AD(40,k);  
   H(k).gain=AD(41,k);  
   H(k).igc=AD(42,k);   
   H(k).igi=AD(43,k);  
   H(k).corr=AD(44,k);   
   H(k).sfs=AD(45,k);   
   H(k).sfe=AD(46,k);   
   H(k).slen=AD(47,k); 
   H(k).styp=AD(48,k);
   H(k).stas=AD(49,k);  
   H(k).stae=AD(50,k);
   H(k).tatyp=AD(51,k);  
   H(k).afilf=AD(52,k); 
   H(k).afils=AD(53,k); 
   H(k).nofilf=AD(54,k);
   H(k).nofils=AD(55,k);
   H(k).lcf=AD(56,k);  
   H(k).hcf=AD(57,k);
   H(k).lcs=AD(58,k);  
   H(k).hcs=AD(59,k);
   H(k).year=AD(60,k);  
   H(k).day= AD(61,k);   
   H(k).hour=AD(62,k); 
   H(k).minute=AD(63,k);
   H(k).sec=AD(64,k);
   H(k).timbas=AD(65,k);
   H(k).trwf=AD(66,k);  
   H(k).grnors=AD(67,k);
   H(k).grnofr=AD(68,k);
   H(k).grnlof=AD(69,k);
   H(k).gaps=AD(70,k);   
   H(k).otrav=AD(71,k);  
   H(k).d1=AD(72,k);     
   H(k).f1=AD(73,k);    
   H(k).d2=AD(74,k);    
   H(k).f2=AD(75,k);    
   H(k).ungpow=AD(76,k);
   H(k).unscale=AD(77,k);
   H(k).ntr=AD(78,k);    
   H(k).mark=AD(79,k); 
%  mp2 = [H.unass]  not used in this prog.        
end 

    size(AD);
    D = AD(79:79+ns,:);

   return

