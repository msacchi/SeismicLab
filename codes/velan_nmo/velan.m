function [S,tau,v] = velan(d,dt,h,vmin,vmax,nv,R,L);
%VELAN: A program to compute velocity spectra.
% 
%  [S,tau,v] = velan(d,dt,h,vmin,vmax,nv,R,L);
%
%  IN   data:      data
%       dt:        sampling interval in secs
%       h:         vector of offsets in meters
%       vmin:      min velocity to scan
%       vmax:      max velocity to scan
%       nv:        number of velocities
%       R:         integer (2,3,4) indicating that the semblance 
%                  is computed every R time samples 
%       L:         the length of the gate of analysis is 2L+1
%
%  OUT  S:         Measure of energy - in this case unnormalized 
%                  cross-correlation in the gate
%       v,tau:     axis vectors to plot the semblance using, for instance, 
%                  imagesc(v,tau,S)
%
%  Reference: Yilmaz, 1987, Seismic Data Processing, SEG Publication
%
%  Example: see va_demo.m
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.



  [nt, nh] = size(d);

  v = linspace(vmin,vmax,nv);

  nv = length(v);
  tau = (0:R:nt-1)*dt;
  ntau = length(tau);
  taper = hamming(2*L+1);
  H = hamming(2*L+1)*ones(nh,1)';

  for it = 1:ntau
  for iv = 1:nv  

  time = sqrt( tau(it)^2 + (h/v(iv)).^2 );

   s = zeros(2*L+1,nh);

    for ig = -L:L;
      ts = time + (ig-1)*dt;

    for ih = 1:nh

   is = ts(ih)/dt+1;
   i1 = floor(is);
   i2 = i1 + 1;

  if i1>=1 & i2<=nt ;
   a = is-i1;
   s(ig+L+1,ih) = (1.-a)*d(i1,ih) + a*d(i2,ih);   % Grab sample with Linear interpolation
  end;

  end
  end

   ss = s.*H;
   s1  = sum( (sum(ss,2)).^2);
   s2  = sum( sum(s.^2));
   S(it,iv) = s1./(s2+0.001);

  end
  end

 S = S/max(max(S));

 return;
    
