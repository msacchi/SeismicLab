function [d]=forward_radon_freq(m,dt,h,p,N,flow,fhigh,ntd);
%FORWARD_RADON_FREQ: Forward linear and parabolic Radon transform.
%                    Freq. domain algorithm
%
%  [d] = forward_radon_freq(m,dt,h,p,N,flow,fhigh,ntd);
%
%  IN   m:     the Radon panel, a matrix m(nt,np)
%       dt:    sampling in sec
%       h(nh): offset or position of traces in mts
%       p(np): ray parameter  to retrieve if N=1
%              curvature of the parabola if N=2
%       N:     N=1 linear tau-p
%              N=2 parabolic tau-p
%       flow:  min freq. in Hz
%       fhigh: max freq. in Hz
%       ntd: number of points of output (Use nt if you don't know)
%
%  OUT  d:     data
%
%  Reference: Hampson, D., 1986, Inverse velocity stacking for multiple elimination,
%             Journal of the CSEG, vol 22, no 1., 44-55.
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.
% 


[nt,nq] = size(m);
nh = length(h);

 if N==2;  
   h=(h/max(abs(h))).^2; 
 end;   

 nfft = 4*2^nextpow2(nt); 

 M = fft(m,nfft,1);
 D = zeros(nfft,nh);
 i = sqrt(-1);

 if1 = floor(flow*dt*nfft)+1;
 if2 = floor(fhigh*dt*nfft)+1;
 
 for ifreq = if1 : if2;
  f = 2.*pi*(ifreq-1)/(nfft*dt);
  L = exp(-i*f*(h'*p));
  x = M(ifreq,:).';
  y = L * x; 
  D(ifreq,:) = y.';
 end

for ifreq = nfft/2+2:nfft 
  D(ifreq,:) = conj(D(nfft-ifreq+2,:));
 end    

d = ifft(D,[],1);

d = d(1:ntd,:); 

return;



