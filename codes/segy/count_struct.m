function [count]=count_struct;
%COUNT_STRUCT: delimiters of header words in the trace header.
%              This is the engine for reading su files
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.

      count.tracl= 0;
      count.tracr= 4;
       count.fldr= 8;
      count.tracf= 12;
         count.ep= 16;
        count.cdp= 20;
       count.cdpt= 24;
       count.trid= 28;
        count.nva= 30;
        count.nhs= 32;
       count.duse= 34;
     count.offset= 36;
      count.gelev= 40;
      count.selev= 44;
     count.sdepth= 48;
       count.gdel= 52;
       count.sdel= 56;
      count.swdep= 60;
      count.gwdep= 64;
     count.scalel= 68;
     count.scalco= 70;
         count.sx= 72;
         count.sy= 76;
         count.gx= 80;
         count.gy= 84;
     count.counit= 88;
      count.wevel= 90;
     count.swevel= 92;
        count.sut= 94;
        count.gut= 96;
      count.sstat= 98;
      count.gstat= 100;
      count.tstat= 102;
       count.laga= 104;
       count.lagb= 106;
      count.delrt= 108;
       count.muts= 110;
       count.mute= 112;
         count.ns= 114;
         count.dt= 116;
       count.gain= 118;
        count.igc= 120;
        count.igi= 122;
       count.corr= 124;
        count.sfs= 126;
        count.sfe= 128;
       count.slen= 130;
       count.styp= 132;
       count.stas= 134;
       count.stae= 136;
      count.tatyp= 138;
      count.afilf= 140;
      count.afils= 142;
     count.nofilf= 144;
     count.nofils= 146;
        count.lcf= 148;
        count.hcf= 150;
        count.lcs= 152;
        count.hcs= 154;
       count.year= 156;
        count.day= 158;
       count.hour= 160;
     count.minute= 162;
        count.sec= 164;
     count.timbas= 166;
       count.trwf= 168;
     count.grnors= 170;
     count.grnofr= 172;
     count.grnlof= 174;
       count.gaps= 176;
      count.otrav= 178;
         count.d1= 180;
         count.f1= 184;
         count.d2= 188;
         count.f2= 192;
     count.ungpow= 196;
    count.unscale= 200;
        count.ntr= 204;
       count.mark= 208;
      count.unass= 210;
      count.trace= 240;

