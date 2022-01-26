%PARABOLIC_MOVEOUT_DEMO: Analysis of residual moveout after NMO.
%
%           -----------------------------------
%              A demo to examine Parabolic 
%            Moveout after  NMO correction
%
%                M.D.Sacchi, SeismicLab
%           -----------------------------------
%
% The residual moveout of a reflection after NMO correction
% can be approximated by a parabolic traveltime curve of the form
% t = tau + q (offset/max(offest))^2, q is the residual moveout
% at far offset. This scripts computes coherence along parabolic
% traveltimes to data that were nmo corrected.
%
% This demo uses the following SeismicLab functions:
% 
%               readsegy.m 
%               nmo.m
%               parabolic_moveout.m
%               seismic.m
 


  clear;
  close all;

% Read a file (su format)

  [D,H] = readsegy('syn_cmp_mult.su');

% Extract offet and sampling interval

  h = [H.offset];
  dtsec = H(1).dt/1000/1000

% Apply constant velocity NMO 

  tnmo = []; 
  vnmo = 3000;
  max_stretch = 40;

  D1 = nmo(D,dtsec,h,tnmo,vnmo,max_stretch);

% Compute measure of coherence as a function
% of parabolic moveout 

  qmin = -0.5;
  qmax = 1.2;
  nq = 120;
  R = 1;
  L = 20;

  [S,tau,q ] = parabolic_moveout(D1,dtsec,h,qmin,qmax,nq,R,L);

% Display results

  figure(1); clf;

  [nt,nx] = size(D);

% Data 

  h1 = figure(1); set(h1,'Position',[250 200 900 1200])

  subplot(131); 
  imagesc(h,[0:1:nt-1]*dtsec,D); 
  xlabel('Offset [m]'); 
  ylabel('Time [s]');
  grid;

% Data after NMO

  subplot(132); 
  imagesc(h,[0:1:nt-1]*dtsec,D1); 
  xlabel('Offset [m]'); 
  grid;

% Energy Spectra versus residual moveout at far offset (q)

  subplot(133); 
  Smax = max(S(:)); 
  imagesc(q,tau,S,[0,0.5*Smax]);  
  xlabel('q (Residual moveout) [s]');
  grid;
  colormap(1-gray);

