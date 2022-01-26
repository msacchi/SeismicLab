%RADON_DEMO_1: Multiple removal with the Parabolic Radon Transform.
%
%
%        -------------------------------------------
%           A demo showing how one can eliminate
%              multiples using the Parabolic 
%                     Radon Transform 
%                   (Synthetic data demo)
%
%                 M.D.Sacchi, SeismicLab
%        -------------------------------------------
%
%
% We read a cmp gather, we apply NMO correction and then
% parabolic Radon demultiple is used to estimate the primaries.
% D = Prim + Mult, an estimate of the primaries is obtained by
% modelling the multiples with the parabolic Radon transform:
% Estimate of Prim = D - Modelled Multiples
%
% This demo uses the following SeismicLab functions:
%
%               readsegy.m
%               nmo.m
%               pradon_demultiple.m
%               inmo.m
%               clip.m
%               seismic.m

  clear;
  close all;


% Read a file (su format). The cmp contains a long period
% multiple of v=1500m/s

  [D,H] = readsegy('syn_cmp_mult.su');

% Extract offet and sampling interval

  h = [H.offset];
  dtsec = H(1).dt/1000/1000;

% Apply constant velocity NMO with the velocity of primaries

  tnmo = [1., 2.]; 
  vnmo = [1500,2000];
  max_stretch = 50;

  D1 = nmo(D,dtsec,h,tnmo,vnmo,max_stretch);

% Use the parabolic Radon Transform to estimate primaries

  qmin = -0.3;
  qmax = .8;
  nq = 60;
  flow = 1.;
  fhigh = 60;
  mu = 1;
  q_cut = 0.01;

  [D1p,M,tau,q] = pradon_demultiple(D1,dtsec,h,qmin,qmax,nq,flow,fhigh,mu,q_cut,'hr');

% Display results with NMO

  h1 = figure(1); set(h1,'Position',[250 200 500 500])
  clipd = 0.9*max(D1(:))*[-1,1];
  clipm = 0.9*max(M(:))*[-1,1];
  subplot(131); imagesc(h,tau,D1,  clipd);  xlabel('Offset [m]'); ylabel('Time [s]');
  subplot(132); imagesc(q,tau,M,   clipm);  xlabel('q (Residual Moveout) [s]')
  subplot(133); imagesc(h,tau,D1p, clipd);  xlabel('Offset [m]')
  colormap(seismic(1));

% Apply Inverse NMO to the model of primaries

  Prim = inmo(D1p,dtsec,h,tnmo,vnmo,max_stretch);

% Display original data and filtered data

  h2 = figure(2); set(h2,'Position',[550 100 500 500])
  clip = 0.9*max(D(:))*[-1,1];
  subplot(121); imagesc(h,tau,D,clip);    xlabel('Offset [m]'); ylabel('Time [s]');
  subplot(122); imagesc(h,tau,Prim,clip); xlabel('Offset [m]');
  colormap(seismic(1));


