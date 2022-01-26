function [segy] = make_empty_header;
%MAKE_EMPTY_HEADER: Make an empty SU-SEGY header.
%
%  [H] = make_empty_header;
%
%  IN   nothing
%
%  OUT  segy: header for one trace with null values
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.

      segy.tracl= 0;
      segy.tracr= 0;
       segy.fldr= 0;
      segy.tracf= 0;
         segy.ep= 0;
        segy.cdp= 0;
       segy.cdpt= 0;
       segy.trid= 0;
        segy.nva= 0;
        segy.nhs= 0;
       segy.duse= 0;
     segy.offset= 0;
      segy.gelev= 0;
      segy.selev= 0;
     segy.sdepth= 0;
       segy.gdel= 0;
       segy.sdel= 0;
      segy.swdep= 0;
      segy.gwdep= 0;
     segy.scalel= 0;
     segy.scalco= 0;
         segy.sx= 0;
         segy.sy= 0;
         segy.gx= 0;
         segy.gy= 0;
     segy.counit= 0;
      segy.wevel= 0;
     segy.swevel= 0;
        segy.sut= 0;
        segy.gut= 0;
      segy.sstat= 0;
      segy.gstat= 0;
      segy.tstat= 0;
       segy.laga= 0;
       segy.lagb= 0;
      segy.delrt= 0;
       segy.muts= 0;
       segy.mute= 0;
         segy.ns= 0;
         segy.dt= 0;
       segy.gain= 0;
        segy.igc= 0;
        segy.igi= 0;
       segy.corr= 0;
        segy.sfs= 0;
        segy.sfe= 0;
       segy.slen= 0;
       segy.styp= 0;
       segy.stas= 0;
       segy.stae= 0;
      segy.tatyp= 0;
      segy.afilf= 0;
      segy.afils= 0;
     segy.nofilf= 0;
     segy.nofils= 0;
        segy.lcf= 0;
        segy.hcf= 0;
        segy.lcs= 0;
        segy.hcs= 0;
       segy.year= 0;
        segy.day= 0;
       segy.hour= 0;
     segy.minute= 0;
        segy.sec= 0;
     segy.timbas= 0;
       segy.trwf= 0;
     segy.grnors= 0;
     segy.grnofr= 0;
     segy.grnlof= 0;
       segy.gaps= 0;
      segy.otrav= 0;
         segy.d1= 0;
         segy.f1= 0;
         segy.d2= 0;
         segy.f2= 0;
     segy.ungpow= 0;
    segy.unscale= 0;
        segy.ntr= 0;
       segy.mark= 0;
      segy.unass= zeros(15,1);
      segy.trace= 0;
