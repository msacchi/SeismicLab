function [Out]=Operator_Radon_Freq(In,Param,flag)
%  [d] = forward_radon_freq(m,dt,h,p,N,flow,fhigh);
%
%  IN   m:     the Radon panel, a matrix m(nt,np)
%       dt:    sampling in sec
%       h(nh): offset or position of traces in mts
%       p(np): ray parameter  to retrieve if N=1
%              curvature of the parabola if N=2
%       N:     N=1 linear tau-p
%              N=2 parabolic tau-p
%       flow:  min freq. in Hz
%       fhigh: max freq. in Hz
%
%  OUT  d:     data
%
%  Reference: Hampson, D., 1986, Inverse velocity stacking for multiple elimination,
%             Journal of the CSEG, vol 22, no 1., 44-55.
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%


   dt = Param.dt;
    h = Param.h;
    p = Param.p;
   nt = Param.nt;
 ntau = Param.ntau;
 nt   = Param.nt;
 flow = Param.flow;
fhigh = Param.fhigh;
N     = Param.N;
Mutes = Param.Mutes; 

nh = length(h);
np = length(p);
nfft = 2*(2^nextpow2(ntau));

if flag == 1;
    m = In;
else
    d = Mutes.*In;
end

if N==2;
    h=(h/max(abs(h))).^2;
end;

ilow  = floor(flow*dt*nfft)+1;
ihigh = floor(fhigh*dt*nfft)+1;
i = sqrt(-1);
A = h'*p;

if flag == 1;
    
    M = fft(m,nfft,1);
    D = zeros(nfft,nh);
    
    for ifreq=ilow:ihigh
        f = 2.*pi*(ifreq-1)/nfft/dt;
        L = exp(i*f*(A));
        x = M(ifreq,:)';
        y = L * x;
        D(ifreq,:) = y';
    end
    
    for k=nfft/2+2:nfft
        D(k,:) = conj(D(nfft-k+2,:));
    end
    
    d = nfft*ifft(D,[],1);
    d = real( d(1:nt,:));
    
    Out = Mutes.*d/(nh*np);
    
else
    
    D = fft(d,nfft,1);
    M = zeros(nfft,np);
    
    for ifreq=ilow:ihigh
        f = 2.*pi*(ifreq-1)/nfft/dt;
        L = exp(i*f*(A));
        y = D(ifreq,:)';
        x = L' * y;
        M(ifreq,:) = x';
    end
    
    for k=nfft/2+2:nfft
        M(k,:) = conj(M(nfft-k+2,:));
    end
    
    m = nfft*ifft(M,[],1);
    m = real( m(1:ntau,:));
    
    
    Out = m/(nh*np);
    
end

return;



