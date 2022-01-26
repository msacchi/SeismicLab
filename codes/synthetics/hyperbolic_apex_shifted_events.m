function [d,h,t] = hyperbolic_apex_shifted_events(dt,f0,tmax,h,h0,tau,v,amp);
%HYPERBOLIC_APEX_SHIFTED_EVENTS: A program to create data containing hyperbolas with shifted apexes 
%
%  [d] = hyperbolic_apex_shifted_events(dt,f0,tmax,h,h0, tau,v,amp,snr,L);
%
%  IN   dt:           sampling interval in secs
%       f0:           central freq. of a Ricker wavelet in Hz
%       tmax:         maximun time of the simulation in secs
%       h:            vector of offsets in meters
%       tau,v,amp,h0: vectors of intercept, rms velocities
%                     ,amplitudes and apexes
%                     (v is in m/s,  tau in s, h0 in m)
%
%  OUT  d:            Data that consist of a superposition of reflections with 
%                     with hyperbolic moveout with apex shifts. 
%       t,h:          time and offset axes 
%
%  Example with default parameters:
%
%    [d,h,t] = hyperbolic_apex_shifted_events; imagesc(h,t,d);
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.



 if nargin == 0
  dt = 2./1000;
  tmax = 1.2;
  h = [-1500:20:1200];
  tau = [0.1,.5,0.8];
  v = [1500,2400,2300];
  h0 = [500,-300,300];
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

% Important: the following lines are needed to have the maximum of the Ricker
% wavelet at the right intercept time

 delay = dt*(floor(nw/2)+1);

 for ifreq=1:nfft/2+1
  w = 2.*pi*(ifreq-1)/nfft/dt;
   for k=1:n_events
    Shift = exp(-i*w*(  sqrt(tau(k)^2 + ((h-h0(k))/v(k)).^2) - delay));
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


