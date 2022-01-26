function [segy] = segy_struct;
%SEGY_STRUCT: Header words for trace SU-SEGY header.
%             This is needed for reading su files
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


      segy.tracl= 'int';
      segy.tracr= 'int';
       segy.fldr= 'int';
      segy.tracf= 'int';
         segy.ep= 'int';
        segy.cdp= 'int';
       segy.cdpt= 'int';
       segy.trid= 'short';
        segy.nva= 'short';
        segy.nhs= 'short';
       segy.duse= 'short';
     segy.offset= 'int';
      segy.gelev= 'int';
      segy.selev= 'int';
     segy.sdepth= 'int';
       segy.gdel= 'int';
       segy.sdel= 'int';
      segy.swdep= 'int';
      segy.gwdep= 'int';
     segy.scalel= 'short';
     segy.scalco= 'short';
         segy.sx= 'int';
         segy.sy= 'int';
         segy.gx= 'int';
         segy.gy= 'int';
     segy.counit= 'short';
      segy.wevel= 'short';
     segy.swevel= 'short';
        segy.sut= 'short';
        segy.gut= 'short';
      segy.sstat= 'short';
      segy.gstat= 'short';
      segy.tstat= 'short';
       segy.laga= 'short';
       segy.lagb= 'short';
      segy.delrt= 'short';
       segy.muts= 'short';
       segy.mute= 'short';
         segy.ns= 'short';
         segy.dt= 'short';
       segy.gain= 'short';
        segy.igc= 'short';
        segy.igi= 'short';
       segy.corr= 'short';
        segy.sfs= 'short';
        segy.sfe= 'short';
       segy.slen= 'short';
       segy.styp= 'short';
       segy.stas= 'short';
       segy.stae= 'short';
      segy.tatyp= 'short';
      segy.afilf= 'short';
      segy.afils= 'short';
     segy.nofilf= 'short';
     segy.nofils= 'short';
        segy.lcf= 'short';
        segy.hcf= 'short';
        segy.lcs= 'short';
        segy.hcs= 'short';
       segy.year= 'short';
        segy.day= 'short';
       segy.hour= 'short';
     segy.minute= 'short';
        segy.sec= 'short';
     segy.timbas= 'short';
       segy.trwf= 'short';
     segy.grnors= 'short';
     segy.grnofr= 'short';
     segy.grnlof= 'short';
       segy.gaps= 'short';
      segy.otrav= 'short';
         segy.d1= 'float';
         segy.f1= 'float';
         segy.d2= 'float';
         segy.f2= 'float';
     segy.ungpow= 'float';
    segy.unscale= 'float';
        segy.ntr= 'int';
       segy.mark= 'short';
      segy.unass= 'short';
      segy.trace= 'float';
