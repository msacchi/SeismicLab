function [d,h,t] = linear_events(dt,f0,tmax,h,tau,p,amp);
%LINEAR_EVENTS: A program to generate t-x data containing linear events.
%
% [d] = linear_events(dt,f0,tmax,h,tau,p,amp);
%
% IN   dt:        sampling interval in secs
%      f0:        central freq. of a Ricker wavelet in Hz
%      tmax:      maximun time of the simulation in secs
%      h:         vector of desire offsets in meters
%      tau,p,amp: vectors of intercept, ray parameter 
%                 and amplitude of each linear event
%                 (p is in sec/m and tau in secs)
%
% OUT  d:         data that consist of a superposition of linear events
%  
% Example:        [d,h,t] = linear_events; imagesc(h,t,d);
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


if nargin == 0
  dt = 4./1000;
  tmax = 1.;
  h = [0:5:5*(55-1)];
  tau = [0.8,0.4,0.33,0.3];
  p = [0.00031,-0.0003,0.0001,0];
  amp = [1.2,-1.,1.,0.10,-1.];
  f0 = 20;
end;
 
 h = abs(h);
 nt = floor(tmax/dt)+1;
 nfft = 4*(2^nextpow2(nt));
 n_events = length(tau);
 nh = length(h);
 wavelet = ricker(f0,dt); 
 nw = length(wavelet);
 W = fft(wavelet,nfft);
 D = zeros(nfft,nh);
 i = sqrt(-1);

 delay = dt*(nw-1)/2;

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

 if nargout>1;
  t = (0:1:nt-1)*dt; 
 end;

 return;
