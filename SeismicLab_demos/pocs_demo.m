%POCS_DEMO: Demo showing POCS reconstruction
%
%
%     -------------------------------------------
%   A demo showing how to run POCS data restoration
%
%      -------------------------------------------
%
%
%               linear_events.m               
%               add_noise.m 
%               pocs.m
%

clear all
close all;


 dt = 0.004;
 f0 = 30;
 nt = 140;
 n1 = 80;
 n2 = 60; 
 dx = [2,2]; 
 amp = [-1, 1    2,    1   , -1];
 p  = [0.04, -0.02;
       0.03,  0.04];
 t0 = [ 0.1, 0.2];
 A  = [-1.0, 1.0];
 snr =9999;
 L = 3; 
% Clearn signal 

 D0 = data_cube([n1,n2],dt,f0, nt,dx,t0,p,A,'parabolic');

% Add noise 

 Dn = add_noise(D0,snr,L);

% Sampling
T = ones(n1,n2);

  for k1 = 1:n1
  for k2 = 1:n2
 p = randn(1,1);
 if p<0.8; 
   Dn(:,k1,k2)=0;
   T(k1,k2)=0.;
end;
   end
  end

 f_low =   0.1;
 f_high = 90.0;
 option = 3;
 perc_i = 99.0;
 perc_f =   .1;
 N = 100; tol = 0.0001;
 
 a = 0.6;
 [Dr,e1,e2,freq] = pocs(Dn,D0,T,dt,f_low,f_high, option, perc_i, perc_f, N, a, tol);


 c = max(abs(D0(:)));
 figure(1); clf;
 set(gcf, 'Position',  [100, 400,1200, 400])
 subplot(131); imagesc(D0(:,:,4),c*[-1,1]);title('Clean signal');
 subplot(132); imagesc(Dn(:,:,4),c*[-1,1]);title('Noisy signal');
 subplot(133); imagesc(Dr(:,:,5),c*[-1,1]);title('Restored signal');
 colormap(seismic(2));

 figure(2); clf;
 set(gcf, 'Position',  [300, 400,1200, 400])
 wigb([squeeze(D0(:,5,:)),squeeze(Dn(:,5,:)),squeeze(Dr(:,5,:))])
 title('Clean - Decimated - Restored')

