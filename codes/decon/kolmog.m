function [w_min] = kolmog(w,type_of_input,L);  
%KOLMOG: Kolmogoroff spectral factorization.
%        Given a wavelet, this function retrieves the minimum 
%        phase wavelet using Kolmogoroff factorization.
%        If the input is a trace, the spectral factorization
%        is applied to the autocorrelation after smoothing.
%
%  [w_min] = kolmog(w,type_of_input,L)
%
%  IN   w:     a wavelet of arbitrary phase if
%              type_of_input = 'w' 
%              or a seismic trace if 
%              type_of_input = 't'
%       L:     lenght of wavelet if type_of_input='t'
%
%  OUT  w_min: a min phase wavelet 
%
%  Reference: Claerbout, 1976, Fundamentals of geophysical data processing 
%     
%  Example:
%
%    w = [1;2;-1;0;0;0;0];
%    wmin = kolmog(w,'w');
%    subplot(221); stem(w);
%    subplot(222); stem(wmin);
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.



 [n1,n2]=size(w);

 if n2~=1; 
  error('Input wavelet or trace must be a column vector');
 end;

if isequal(type_of_input,'w')

 nw = length(w); 
 nfft = 8*( 2^nextpow2(nw));
 W = log ( abs(fft(w,nfft)) +0.00001);
 W = ifft(W);
 for k=nfft/2+2:nfft; 
  W(k)=0.;
 end;
 W = 2.*W;
 W(1) =W(1)/2.;
 W = exp(fft(W)) ;
 w_min = real(ifft(W));
 w_min = w_min(1:nw); 

else;

 nt = length(w);   
 nfft = 8*( 2^nextpow2(nt));
 nw = L;
 A = xcorr(w,w,L);
 A = A.*hamming(2*L+1);
 W = log ( sqrt(abs(fft(A,nfft))) +0.00001);
 W = ifft(W);
 for k=nfft/2+2:nfft; W(k,1)=0.;end;
 W = 2.*W;
 W(1,1) =W(1,1)/2.;
 W = exp(fft(W)) ;
 w_min = real(ifft(W));
 w_min = w_min(1:nw,1); 
 w_min_max = max(abs(w_min));
 w_min = w_min/w_min_max;
end;

return

