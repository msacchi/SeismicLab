function [P,f,w,tw]  = smooth_spectrum(d,dt,L,io);
%SMOOTH_SPECTRUM: Power spectrum estimation by smoothing the periodogram.
%                 For more than one trace provides the average spectrum
%                 followed by smoothing.
%
%  [P,f,w,tw]  = smooth_spectrum(d,dt,L,io);
%
%  IN     d:  data  array  data(t,x) or data(t,x,y) or data(t,x,y,h) etc 
%         dt: sampling interval in secs
%         L:  Lenght of the freq. smoothing operator
%             L=0  means no snoothing 
%         io: 'db' in db, 'li' for linear scale 
%
%  OUT    P:  normalized smoothed power spectrum in linear or dB scale.
%         f:  frequency axis 
%         w:  wavelet (zero phase with Power spectrum P)
%         tw: time axis to plot wavelet
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


 Ndims = ndims(d); 
 N     = size(d);
 nt    = N(1);

% reshape any cube into a matrix 

 aux  = reshape(d,N(1),prod(N(2:end)));

 wind = hamming(2*L+1);

 nf = max(2*2^nextpow2(nt),2048);

 f = 1:1:nf/2+1;
 f = (f-1)/dt/nf;      % Freq. axis in Hz 
 D = fft(aux,nf,1);
 D = sum(abs(D).^2,2);
 D = conv(D,wind);     % Smooth 
 N = length(D);
 D = D(L+1:N-L);
 A = sqrt(D);
 D = D(1:nf/2+1,1);
 D = D/max(D);

if(io=='db'); 
 P = 10*log10(D);      % Power spectrum in dB
 I = find(P<-40); P(I)=-40;
else
P = D;
end

if (nargout>2); 

% Estimate a zero phase wavelet    
% For the length I am assuming a wavelet with central
% freq. 30Hz

  f0 = 30;

  L = 3/(f0*dt);
  w = real(fftshift( ifft(A)));
  w = w(nf/2+1-L:nf/2+1+L,1);
  w = w.*hamming(2*L+1);
  w = w/max(abs(w));
  tw = [-L:1:L]*dt;

end;

