function [d,h,t] = hyperbolic_events(dt,f0,tmax,h,tau,v,amp);
%HYPERBOLIC_EVENTS: A program to generate data containing hyperbolas.
%
%  [d] = hyperbolic_events(dt,f0,tmax,h,tau,v,amp);
%
%  IN   dt:        sampling interval in secs
%       f0:        central freq. of a Ricker wavelet in Hz
%       tmax:      maximun time of the simulation in secs
%       h:         vector of offsets in meters
%       tau,v,amp: vectors of intercept, rms velocities
%                  and amplitude of each linear event
%                  (v is in m/s and tau in secs)
%
%  OUT  d:         Data that consist of a superposition of reflections
%                  with hyerbolic  moveout (no avo)
%       t,h:       time and offset axes 
%
%  Example with default parameters:
%
%    [d,h,t] = hyperbolic_events; imagesc(h,t,d);
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.



 if nargin == 0
  dt = 4./1000;
  tmax = 1.2;
  h = [20:20:1000];
  tau = [0.3,.5,0.8];
  v = [2000,2500,3000];
  amp = [1., -1.,1];
  f0 = 20;
 end;
 
 nt = floor(tmax/dt)+1;
 nfft = 4*(2^nextpow2(nt));
 n_events = length(tau);
 nh = length(h);
 wavelet = ricker(f0,dt); 
 nw = length(wavelet);
 W = fft(wavelet,nfft);
 D = zeros(nfft,nh);
 i = sqrt(-1);

% Important: the following lines area needed to have the maximum of the Ricker
% wavelet at the right intercept time

 delay = dt*(floor(nw/2)+1);

 for ifreq=1:nfft/2+1
  w = 2.*pi*(ifreq-1)/nfft/dt;
   for k=1:n_events
    Shift = exp(-i*w*(  sqrt(tau(k)^2 + (h/v(k)).^2) - delay));
   D(ifreq,:) = D(ifreq,:) +amp(k)* W(ifreq)*Shift;
  end
 end

% Apply w-domain symmetries

 for ifreq=2:nfft/2
  D(nfft+2-ifreq,:) = conj(D(ifreq,:));
 end 

 d = ifft(D,[],1);
 d = real(d(1:nt,:));

 if nargout>1;
  t = (0:1:nt-1)*dt;
 end;

 return;
