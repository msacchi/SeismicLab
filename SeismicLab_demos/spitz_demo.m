%SPITZ_DEMO: Example showing Spitz fx interpolation
%
%        -----------------------------------
%         A demo for Spitz FX interpolation
%
%            M.Naghizadeh and M.D.Sacchi
%                    SeismicLab
%        -----------------------------------
% 
% Spitz FX interpolation is designed for interpolating
% events with linear moveout. It can tolarate, however,
% curvature as shown in this example. In general,
% the method should be applied in small overlapping
% windows.
%
% This demo uses the following SeismicLab functions:
%
%               readsegy.m 
%               spitz_fx_interpolation.m
%               wigb.m
%

clear;
close all;
clc

% Read a file (su format)

  [d,H] = readsegy('gom_cdp_nmo.su');

% Extract offet and sampling interval 

  h = [H.offset];
  dt = H(1).dt/1000/1000;

% Use a subset of the data

  n0 = 400;
  d = d(n0:end,:);

% Time axis for figures

  [nt,nh] = size(d);
  t0 = (n0-1)*dt;
  taxis = t0+[0:1:nt-1]*dt;
  
% New offset coordinates after interpolation

  h_interp = linspace(max(h),min(h),2*nh-1);

% Parameters for interpolation

  npf = 25;            % Length of the pef
  pre1 = 1.0;          % percentage of pre-whitening for estimating the pef
  pre2 = 1.0;          % percentage of pre-whitening for estimating the data
  flow = 0.1;          % first frequency for the f-x process (Hz)
  fhigh = 90;          % last frequncy (Hz)

  [di] = spitz_fx_interpolation(d,dt,1,pre1,pre2,flow,fhigh);
  
% Figures - first with wiggles

  h1 = figure(1); set(h1,'Position',[100 100 800 900])
  subplot(121); wigb(d,2,h,taxis); title('Orignal');
  subplot(122); wigb(di,2,h_interp,taxis); title('After Spitz FX Interpolation');

% Figures - images

  h2 = figure(2); set(h2,'Position',[300 300 800 900])
  subplot(121); imagesc(h,taxis,clip(d,60,60)); title('Orignal');
  subplot(122); imagesc(h_interp,taxis,clip(di,60,60)); title('After Spitz FX Interpolation');
  colormap(seismic);
 

  
