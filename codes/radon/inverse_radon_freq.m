function [m] = inverse_radon_freq(d,dt,h,q,N,flow,fhigh,mu,sol);
%INVERSE_RADON_FREQ: Inverse linear or parabolic Radon transform. 
%                    Frequency domain alg.
%
%  [m] = inverse_radon_freq(d,dt,h,q,N,flow,fhigh,mu,sol)
% 
%  IN   d:     seismic traces   
%       dt:    sampling in sec
%       h(nh): offset or position of traces in meters
%       q(nq): ray parameters  if N=1
%              residual moveout at far offset if N=2
%       N:     N=1 Linear tau-p  
%              N=2 Parabolic tau-p
%       flow:  freq.  where the inversion starts in HZ (> 0Hz)
%       fhigh: freq.  where the inversion ends in HZ (< Nyquist) 
%       mu:    regularization parameter 
%       sol:   'ls' is Hampson (1986) damped least-squares solution
%              'hr' is Sacchi & Ulrych (1995) high-res (Sparse) solution
%
%  OUT  m:     the linear or parabolic tau-p panel
%
%  Reference: Hampson, D., 1986, Inverse velocity stacking for multiple elimination
%             Journal of the CSEG, vol 22, no 1., 44-55.
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


 [nt,nh] = size(d);
 nq = length(q);

 if N==2; h=h/max(abs(h));end;  

 nfft = 4*2^(nextpow2(nt)+1);
 if1 = floor(flow*dt*nfft)+1; 
 if2 = floor(fhigh*dt*nfft)+1; 

 D = fft(d,nfft, 1); 
 M = zeros(nfft,nq);

 i = sqrt(-1);

% Minimum norm-style solution 

 for ifreq= if1 : if2;
  w = 2.*pi*(ifreq-1)/(nfft*dt);
  L = exp(-i*w*(h.^N)'*q); 
  y = D(ifreq,:).';
   if isequal(sol,'hr'); x = solver_hr_radon(L,y,mu); end 
   if isequal(sol,'ls'); x = solver_ls_radon(L,y,mu); end 
  M(ifreq,:) = x.';
 end

 for ifreq = nfft/2+2:nfft
  M(ifreq,:) = conj(M(nfft-ifreq+2,:));
 end

 m = ifft(M,[],1);
 m = m(1:nt,:);

return


function x = solver_ls_radon(L,y,mu);
% Least-squares solver for the freq. domain Radon transform.
% For stability I am adding damping via the tradeoff parameter mu

[nh,np]=size(L); 

% Regularized least-squares solution with damping. 

if nh<=np; 
  A = L*L' + mu*eye(nh);
  u = A\y;
  x  =  L'*u;
else;
  x = (L'*L + mu*eye(np))\(L'*y);
end

return

function x = solver_hr_radon(L,y,mu);
% Least-squares solver for the freq. domain Radon transform.
% For stability I am adding Cauchy like regularization with tradeoff  parameter mu
% This is often called the high-resolution solution or sparse solution of 
% the frequency domain Radon Transform. 

[nh,np]=size(L); 

% Regularized least-squares solution with sparsity promoting term
% This is the basic IRLS code in Sacchi & Ulrych 1995. 

if nh<=np; 

 Q = eye(np); 
 for k = 1:3
  A = L*Q*L' + mu*eye(nh);
  u = A\y;
   x  =  Q*L'*u;
   q = (abs(x).^2 + 0.001); 
   Q = diag(q); 
 end

else;

 Q = eye(np); 
 for k = 1:3
   x = (L'*L + mu*Q)\(L'*y);
   q = 1./(abs(x).^2 + 0.001); 
   Q = diag(q); 
 end

end

return
