function [S,k,f] = fk_spectra(d,dt,dx,L);
%FK_SPECTRA: FK spectrum of a seismic gather.
%
%  [S,k,f] = fk_spectra(d,dt,dx,L);
%
%  IN   d:      data (traces in columns) 
%       dt:     time interval
%       dx:     spatial increment between traces 
%       L:      apply spectral smoothing using a separable
%               2D Hamming window of LxL samples
%
%  OUT  S:      FK spectrum
%       f:      freq axis in Hz
%       k:      wave-number spectrum in cylces/m (if dx is in meters)
%
%  Note: when plotting spectra (S)  use log(S) or S.^alpha (alpha=0.1-0.3) to
%        increase the visibility of small events 
%
%  Example: 
%
%    [d,h,t] = linear_events; dt = t(2)-t(1); dx = h(2)-h(1); 
%    [S,k,f] = fk_spectra(d,dt,dx,6); imagesc(k,f,S);
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License. 
    


[nt,nx]=size(d);

 nk = 4*(2^nextpow2(nx));
 nf = 4*(2^nextpow2(nt));

 S = fftshift( abs(fft2(d,nf,nk)) );
 S = S.^2;
 H = hamming(L)*hamming(L)';
 H = H/sum(H(:));
 S = conv2(S,H,'same');
 S = S(nf/2:nf,:);
 
 f = (0:1:nf/2)/nf/dt;
 k = (-nk/2+1:1:nk/2)/nk/dx;

 return;
