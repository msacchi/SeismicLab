function [w, r,d0,d] = make_section(nt,nx,nJ,SNR,seed1,seed2,rho);
% make a synthetic data set
%
% nt,nx:  size of the stack section
% SNR:    desired S/R ratio
% seed1,seed2: seeds for random number generators
%
% w:  wavelet 
% r:  reflectivity
% d0: data (no noise)
% d:  data corrupted by noise 
%

 dt = 2/1000;               % Samplig interval
 f0 = 40;                   % Central frequency of wavelet
 c = 50;                    % Wavelet constant phase rotation 

 randn('state',seed1);
 rand('state',seed2);

 w = phase_correction( ricker(f0,dt), c);

 r = zeros(nt,nx); 

 J = 10 + floor(rand(nJ,1)*(nt-10))
 a = rand(nJ,1)*2-1; 

 r(J,1)= a;
 for k = 2:nx
 pert = round(rand(nJ,1)*2-1);
 J = J + pert;
 r(J,k) = a+0.06*randn(size(a));
 end;

 d0 = conv2(r,w,'same');

 sed = std(d0(:));
 noise = randn(size(d0));
 noise = conv2(noise,hamming(3),'same'); 
 sen  = std(noise(:));

 num = (sed/sen)/SNR

max(d0(:))
 noise = noise*num;


 d = d0 + noise;
 

return;

