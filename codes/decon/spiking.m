function [f,o] = spiking(d,NF,mu);  
%SPIKING: Spiking deconvolution using Levinson's recursion.
%         Statitical deconvolution of seismic traces assuming
%         white reflectivty and minimun phase source wavelet. 
%
%  [f,o] = spiking(d,NF,mu)
%
%  IN   d: data (traces are columns)
%      NF: lenght of the spiking operator
%      mu: prewhitening in percentage  
%
%  OUT  f: the filter
%       o: the ouput or convolution of the data with 
%          the filter (adjusted to the length of the
%          input data and normalized).
%
%  Note: We assume a minimum phase wavelet, we also assume
%        that the reflectivity is a white process. The latter
%        allows one to estimate the autocorrelation of
%        the wavelet from the autocorrelation of the trace.
%
%  Reference: Robinson and Treitel, 1980, Geophysical Signal Analysis, Prentice Hall
%
%  Note: some clarity was lost in order to use Matlab function "levinson" 
%
%  Example:
%      Make a reflectivity and minimum phase wavelet, then
%      convolve them to obtain the seismic trace. 
%      Then, use spiking decon to recover the reflectivity.
%
%      r = bernoulli_refl(200,[0.9,1]); 
%      w = conv([2,-1]',hamming(7)); w = kolmog(w,'w'); 
%      s = conv(w,r); 
%      [f,r_est] = spiking(s,10,0.1);
%      wigb([s,r_est]);
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


 NF = NF - 1;
 [ns,ntraces] = size(d);
 dmax = max(max(abs(d)));

 R = xcorr(d(:,1),d(:,1),NF);       % Compute data autocorrelation for trace 1

 if ntraces>1;
  Ra = R;
   for k=2:ntraces;
    R = xcorr(d(:,k),d(:,k),NF);    % Compute data mean autocorrelation
    Ra = Ra + R;
   end;
  R = Ra/ntraces;
 end;

 Rs = R(:,1).*hamming(NF*2+1);
 r = Rs(NF+1:2*NF+1,1);
 r(1,1) = r(1,1)*(1 + mu/100.);     % Add pre-whitening for stability

 [f] = levinson(r,NF);              % Fast inversion of Toeplitz system 
                                    % via Levinson's recursive algorithm

 f = [f'];                          % I like column vectors

 if nargout == 2
  o = conv2(d,f);          
  o = o(1:ns,:);
  omax = max(max(abs(o)));
  o = o * dmax/omax;
 end

return

