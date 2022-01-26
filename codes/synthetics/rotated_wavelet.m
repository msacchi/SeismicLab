function [w,t,A,P,freq] = rotated_wavelet(dt,fl,fh,c);
%ROTATED_WAVELET: Band-limited wavelet with phase rotation c.
%
%  [w,t,A,P,freq] = rotated_wavelet(dt,fl,fh,c);
%
%  IN     dt:   sampling interval in sec
%         fl:   min freq. in Hz
%         fh:   max freq. in Hz
%         c:    rotation in degs
%
%  OUT    w:    wavelet (column)
%         t:    time axis in secs
%         A:    Amplitude spectrum
%         P:    phase spectrum in degs
%         freq: freq. axis in Hz
%
%  This is an analytical construction with a box-car function. For a wavelet
%  with trapeziodal amplitude spectrum see trapezoidal_wavelet.m
% 
%     
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.



% Define expected length of the wavelet 

 fc = (fh-fl)/2
 L = 4*floor(2.5/fc/dt);

 t = dt*[-L:1:L]';
 B = 2*pi*fh;
 b = 2*pi*fl;
 c = c*pi/180;

 w = sin(c+B*t)-sin(c+b*t);

% avoid 0/0 (Use L'Hopital rule)

 I = find(t==0);
 t(I)=99999;

 w = w./(t*pi);
 w(I) = (B/pi)*cos(c)-(b/pi)*cos(c);
 t(I) = 0.;

% Normalize with dt to get unit amplitude spectrum 
% and smooth with a Hamming window

 w = dt*w.*hamming(length(w));

 nw = 2*L+1;  nh=L+1;

 nf = 4*2^nextpow2(nw);

% Take into account that the wavelet in non-causual
% The followin will remove linear phase shift and
% show the actual phase of the wavelet

 ww = [w(nh:nw,1);zeros(nf-nw,1);w(1:nh-1,1)];

 W = fft(ww);
 M = length(W)/2+1;
 A = abs(W(1:M,1));
 P = (180/pi)*(angle(W(1:M,1)));
 Kh = floor(fh*nf*dt)+1;
 Kl = floor(fl*nf*dt)+1;
 P(Kh+1:M,1)=0;
 P(1:Kl,1)=0;

% freq axis in Hz

 freq = (0:1:M-1)/dt/nf;


