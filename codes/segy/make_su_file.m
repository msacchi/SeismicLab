function  [H] = make_su_file(filename,D,dtsec,offset);
%MAKE_SU_FILE: Write an su file to filename.
%
%  make_su_file(filename,D,dtsec,offset);
%
%  IN   filename: name of file where the data is saved
%       D:        data (traces are columns)
%       dtsec:    dt in secs
%       offset:   offset (integers)
% 
%  OUT  H:        header in structure 
%
%  This program only assigns very basic header words
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


[nt,nx]=size(D);

for k=1:nx,

 H(k)= make_empty_header;
 H(k).dt=floor(dtsec*1000*1000);
 H(k).d1=dtsec;
 H(k).ns=nt;
 H(k).trace =  D(:,k);
 H(k).tracl = k;
 H(k).offset = offset(k);

end;

writesegy(filename,D,H);
