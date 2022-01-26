function [prim,m,tau,q] = pradon_demultiple(d,dt,h,qmin,qmax,nq,flow,fhigh,mu,q_cut,sol);
%PRADON_DEMULTIPLE: Multiple removal using parabolic Radon Transform.
%
%  [prim,m,tau,q] = pradon_demultiple(d,dt,h,qmin,qmax,nq,flow,fhigh,mu,q_cut,sol);
%
%  IN   d:     data   (d(nt,nh))
%       dt:    sampling interval in secs
%       h:     offset in meters  (h(nh))
%       qmin:  min residual moveout at far offset (secs)
%       qmax:  max residual moveout at far offset (secs)
%       nq:    number of samples of the residual moveout axis
%       flow:  min freq to process (Hz)
%       flow:  max freq to process (Hz)
%       mu:    trade-off parameter for LS solution used
%              to retrieve the Radon panel
%       q_cut: keep contributions from q_cut to qmax,
%              this is to estimate the multiples that 
%              are removed from the primaries
%       sol:   'ls' is Hampson (1986) damped least-squares solution 
%              'hr' is Sacchi & Ulrych (1995) high-res (Sparse) solution
%
%  OUT  prim: primaries obtained by removing from the data the
%             multiples modelled  with the Radon transform
%       m:    panel with the Radon transform  (m(nt,nq))
%       tau:  vertical axis for the Radon panel (tau(nt))
%       q:    horizontal axis for the Radon panel (q(nq))
%
%  References: - Hampson, D., 1986, Inverse velocity stacking for multiple 
%                elimination, Journal of the CSEG, 22(1), 44-55.
%              - Sacchi M.D. and T.J. Ulrych, 1995, High-resolution velocity gathers and 
%                offset space reconstruction, Geophysics, 60(4), 939-1278
%
%  Example: see radon_demo.m 
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.



  [nt,nx] = size(d);
  dq = (qmax-qmin)/(nq-1);
  q = qmin+dq*[0:1:nq-1];

  N = 2;                % parabolic tranform

% Transform from t-x to tau-q

  [m] = inverse_radon_freq(d,dt,h,q,N,flow,fhigh,mu,sol);

  nq = length(q);
  iq_cut = floor((q_cut-qmin)/dq)+1;
  mc = m;

  mc(:,1:iq_cut) = 0;   % Keep multiples in the Radon panel

% Transform from tau-q to t-x

  [dm] = forward_radon_freq(mc,dt,h,q,N,flow,fhigh,nt);

  prim = d-dm;

  tau = (0:1:nt-1)*dt;

  return;
