%VA_DEMO: Test velocity analysis program.
%
%
%          -----------------------------------
%           A demo for VA (Velocity Analysis)
%
%                 M.D.Sacchi, SeismicLab
%           -----------------------------------
%
% This script shows how one can compute the un-normalized
% cross-correlation (a measure of coherency) along hyperbolic paths 
% to extract stacking velocities from cmp/cdp gathers.
%
% This demo uses the following SeismicLab functions:
%
%                readsegy.m
%                velan.m
%

  clear;
  close all;

% Read a file (su format)

  [D,H] = readsegy('syn_cmp.su');

% Extract offset and sampling interval

  h = [H.offset];
  dtsec = H(1).dt/1000/1000

% Define parameters for VA
  
  vmin = 1000;
  vmax = 4000;
  nv = 150;
  R = 1;
  L = 15;

% Compute the measure of coherence S 
% (un-normalized cross-correlation)

  [S,tau,v] = velan(D,dtsec,h,vmin,vmax,nv,R,L);

% Display directly with MATLAB

  figure(1); clf;

  [nt,nx] = size(D);

  subplot(121); 
  imagesc(h,[0:1:nt-1]*dtsec,D); 
  xlabel('Offset [m]'); 
  ylabel('Time [s]');

  subplot(122); 
  imagesc(v,tau,S,[-0.1,0.9]);  
  grid; 
  xlabel('V [m/s]');

  colormap(seismic);

