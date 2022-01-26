function [S,tau,q] = parabolic_moveout(d,dt,h,qmin,qmax,nq,R,L);
%PARABOLIC_MOVEOUT: A program to compute parabolic spectra.
%                   This is a display of energy versus residual moveout
%                   at far offset. Input data is a  NMO-ed corrected
%                   gather.
% 
%  [S,tau,q] = parabolic_moveout(d,dt,h,qmin,qmax,nq,R,L);
%
%  IN   data:      data
%       dt:        sampling interval in secs
%       h:         vector of offsets in meters
%       qmin:      min residual moveout at far offset trace (secs)
%       qmax:      max residual moveout at far offset trace (secs)
%       nq:        number of velocities
%       R:         integer (2,3,4) indicating that the samblance 
%                  is computed every R time samples (use for efficiency)
%       L:         length of the semblance gate is 2L+1
%
%  OUT  S:         Measure of energy - in this case unnormalized cross-correlation in the gate
%       q,tau:     axes to plot the semblance using, for instance,  imagesc(q,tau,S)
%
%
%  Example: see parabolic_moveout_demo.m
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.




  [nt, nh] = size(d);

  q = linspace(qmin,qmax,nq);

  nq = length(q);
  tau = (0:R:nt-1)*dt;
  ntau = length(tau);
  taper = hamming(2*L+1);
  H = hamming(2*L+1)*ones(nh,1)';

  hmax = max(abs(h));

  for it = 1:ntau
  for iq = 1:nq  

  time =  tau(it) + q(iq)*(h/hmax).^2 ;

  s = zeros(2*L+1,nh);

  for ig = -L:L;
  ts = time + (ig-1)*dt;

  for ih = 1:nh

   is = ts(ih)/dt+1;
   i1 = floor(is);
   i2 = i1 + 1;

  if i1>=1 & i2<=nt ;
   a = is-i1;
   s(ig+L+1,ih) =  (1.-a)*d(i1,ih) + a*d(i2,ih) ;
  end;

  end
  end

   s = s.*H;

   s1  = sum( (sum(s,2)).^2);
   s2  = sum( sum(s.^2));
   S(it,iq) = s1/(s2+0.00000001);


  end
  end


 S = S/max(max(S));
 return;
    
