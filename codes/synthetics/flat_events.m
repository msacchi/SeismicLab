function [d,h,t] = flat_events(snr);
%FLAT_EVENTS: A program to generate data containing flat events.
%
%  [d] = flat_events(snr)
%
%  IN   snr:       signal-to-noise ration
%
%  OUT  d:         superposition flat events + noise
%  
%  Example:     
%
%    [d,h,t] = flat_events(snr); imagesc(h,t,d);
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


  dt = 2./1000;
  tmax = 0.8;
  h = [0:10:10*(90-1)];
  tau = [0.1,0.2,0.3,0.6];
  p = [0.0,-0.0,0.,-0.];
  amp = [1.2,-1.,1.,1];
  f0 = 20;
  L = 5;
 
 nt = floor(tmax/dt)+1;
 nfft = 4*(2^nextpow2(nt));
 n_events = length(tau);
 nh = length(h);
 wavelet = ricker(f0,dt); 
 nw = length(wavelet);
 W = fft(wavelet,nfft);
 D = zeros(nfft,nh);
 i = sqrt(-1);

 delay = dt*(floor(nw/2)+1);

 for ifreq=1:nfft/2+1
  w = 2.*pi*(ifreq-1)/nfft/dt;
   for k=1:n_events
    Shift = exp(-i*w*(tau(k)+h*p(k)-delay));
   D(ifreq,:) = D(ifreq,:) + amp(k)* W(ifreq)*Shift;
  end
 end

% w-domain symmetries

 for ifreq=2:nfft/2
  D(nfft+2-ifreq,:) = conj(D(ifreq,:));
 end 

 d = ifft(D,[],1);
 d = real(d(1:nt,:));

% My definition of snr = (Max Amp of Clean Data)/(Max Amp of Noise)

 dmax  = max(max(abs(d)));
 op = hamming(L);
 Noise = conv2(randn(size(d)),op,'same');

 Noisemax = max(max(abs(Noise)));

 d = d + Noise*(dmax/Noisemax)/snr;

 if nargout>1;
  t = (0:1:nt-1)*dt; 
 end;

 return;
