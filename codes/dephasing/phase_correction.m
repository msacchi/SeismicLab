function dout = phase_correction(din,c);
%PHASE_CORRECTION: Apply a constant phase correction to seismic data.
%              
%  [dout] = phase_correction(din,c);
%
%  IN   din:      Seismic data (traces in columns);
%       c:        Constant Phase rotation in degrees
%
%  OUT  dout:     data after correction
%
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


  c = c*pi/180;
  [nt,nx]=size(din);
  nf = 2^nextpow2(nt);
  [nt,nx]=size(din);

  Din = fft(din,nf,1);

  Phase = exp(-i*[0;-c*ones(nf/2-1,1);0;c*ones(nf/2-1,1)]);

   for k=1:nx;
     Din(:,k) = Din(:,k).*Phase;
   end

 dout = ifft(Din,nf,1);
 dout = real(dout(1:nt,:));

 return;
