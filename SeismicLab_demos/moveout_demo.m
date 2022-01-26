%MOVEOUT_DEMO: NMO and Inverse NMO demo.
%
%
%      -----------------------------------
%        A demo for NMO and Inverse NMO
%
%           M.D.Sacchi, SeismicLab
%      -----------------------------------
%
% Show nmo and inmo capabilties and how one can create
% a su file from matlab and display it using SU (suxwigb)
% from Matlab.
%
% This demo uses the following SeismicLab functions:
%
%               readsegy.m 
%               nmo.m
%               inmo.m
%               writesegy.m
%               seismic.m
%               suxwigb from SU
%

clear all
close all;


% Read a file (su format)

  [D,H] = readsegy('syn_cmp.su');

% Extract offet and sampling interval

  h = [H.offset];
  dtsec = H(1).dt/1000/1000

% Define t-v pais for NMO

  tnmo = [0.5,  1.22, 1.65];
  vnmo = [2000, 5000, 2500];

  Max_Stretch = 400;

% NMO

  D1 = nmo(D,dtsec,h,tnmo,vnmo,Max_Stretch);

% Inverse NMO

  D2 = inmo(D1,dtsec,h,tnmo,vnmo,Max_Stretch);

% Write results in su format

  writesegy('syn_cmp_nmo.su', D1,H);
  writesegy('syn_cmp_inmo.su',D2,H);

% Display results with SU (remove the next 2 lines 
% if SU is not installed in your system)

 !suxwigb<syn_cmp_nmo.su perc=99 key=offset xbox=200&
 !suxwigb<syn_cmp_inmo.su perc=99 key=offset xbox=400&

% Alternatevely, one can display directly with MATLAB

  figure(1); clf;

  [nt,nx] = size(D);
  imagesc([1:1:2*nx],[0:1:nt-1]*dtsec,[D,D1]);
  colormap(seismic);
