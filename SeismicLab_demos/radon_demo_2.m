%RADON_DEMO_2: Multiple removal with the Parabolic Radon Transform.
%
%        -------------------------------------------
%           A demo showing how one can eliminate
%              multiples using the Parabolic 
%                     Radon Transform 
%                    (real data demo)
%
%                 M.D.Sacchi, SeismicLab
%        -------------------------------------------
%
%
% Parabolic Radon demultiple is used to estimate the primaries.
% D = Prim + Mult, an estimate of the primaries is obtained by
% modelling the multiples with the parabolic Radon transform:
% Estimate of Prim = D - Modelled Multiples
%
% This demo uses the following SeismicLab functions:
%
%               readsegy.m
%               pradon_demultiple.m
%               clip.m
%               seismic.m

  clear;
  clc;
  close all;
  fprintf(' ----------------------------------------- \n')
  fprintf(' Demo Parabolic Radon Transform Demultiple \n')
  fprintf(' ----------------------------------------- \n')


% Read data - a nmo-ed corrected marine gather
% from the Gulf of Mexico

  [D,H] = readsegy('gom_cdp_nmo.su');

% We will use a subset of the data

  D = D(600:1200,:);

% Extract areas that were muted

  Mutes = ones(size(D)); 
  I = find(D==0); 
  Mutes(I)=0;

% Extract offet and sampling interval

   h = [H.offset];
  dt = H(1).dt/1000/1000;

% Define parameters for Radon Demultiple

  qmin = -0.9;
  qmax =  1.2;
  nq = 180;
  flow = 0.1;
  fhigh = 90.0;
  mu = 10.2;
  q_cut = .05;



  [P,M,tau,q] = pradon_demultiple(D,dt,h,qmin,qmax,nq,flow,fhigh,mu,q_cut,'ls');

% Apply data mutes back to estimate of primaries

  P = P.*Mutes;

% Display results

  h1 = figure(1);
  set(h1,'Position',[10 10 500 500])
  t0 = (600-1)*dt; 
  subplot(131); pimage(h,t0+tau,clip(D,50));
  xlabel('Offset [ft]');
  ylabel('Time [s]');
  subplot(132); pimage(q,t0+tau,clip(M,50));
  xlabel('Residual moveout [s]');
  subplot(133); pimage(h,t0+tau,clip(P,50));
  xlabel('Offset [ft]');
  colormap(seismic(1));

  [Sd,tau,q] = parabolic_moveout(D, dt, h, qmin, qmax, nq, 2, 15);
  [Sp,tau,q] = parabolic_moveout(P, dt, h, qmin, qmax, nq, 2, 15);

  h2 = figure(2);
  set(h2,'Position',[50 10 500 500])
  subplot(121); pimage(q,t0+tau,Sd);
  ylabel('Time [s]');
  xlabel('Residual moveout [s]');
  subplot(122); pimage(q,t0+tau,Sp);
  ylabel('Time [s]');
  xlabel('Residual moveout [s]');
  colormap(seismic(1));


