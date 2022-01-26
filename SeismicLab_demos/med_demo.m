%MED_DEMO: Demo for minimum entropy deconvolution (MED)
%
% This demo uses the following function
% 
%               hyperbolic_events.m
%               add_noise.m 
%               pocs.m
%               med.m
%


clear
clc
close all

  dt = 4./1000;
  tmax = 1.2;
  h = [0:12:1000];
  tau =      [0.1, 0.22   0.24,   0.39, 0.44,  0.59, 0.8,  0.82];
  v   = 1000*[1.9, 1.20   1.19,   2.39, 2.44,  2.59, 2.8,  2.09];
  amp =      [1.1, 1.00  -1.01,  -1.09, 0.44, -0.59, 0.8, -1.09];
  f0 = 20;


[s,h]=readsegy('gom_near_offsets.su'); 

 %[s] = hyperbolic_events(dt,f0,tmax,h,tau,v,amp);
  %s = add_noise(s,10,3);

wbp = [1]';

dt = 0.004
Nf = 21;
mu = 0.01;
Updates = 24;



 s = s(400:800,:);
 [filter,tf,x,tx,Med_Norm] = med(wbp,s,dt,Nf,mu,Updates,@funalpha,@dfunalpha,3);

 [Ps,f] = smooth_spectrum(s,dt,5,'li');
 [Px,f] = smooth_spectrum(x,dt,5,'li');
 


% Display results with NMO

  nx = size(x,2);
  h1 = figure(1); set(h1,'Position',[250 200 900 500])
  clipd = 0.2*max(s(:))*[-1,1];
  subplot(121); imagesc([1:1:nx],tx,s,  clipd);  xlabel('Offset [m]'); ylabel('Time [s]');
  subplot(122); imagesc([1:1:nx],tx,x,  clipd);  xlabel('Offset [m]'); ylabel('Time [s]');
  colormap(seismic(1));

 

  h2 = figure(2); set(h2,'Position',[450 200 400 500])
  clipd = 0.2*max(s(:))*[-1,1];
  plot(f,Ps,f,Px);

 [w,o] = ls_inv_filter(filter,40,20,0.001)

 

