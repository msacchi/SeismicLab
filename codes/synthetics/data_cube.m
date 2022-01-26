function  [d] = data_cube(N,dt,f0,nt,dx,t0,p,A,opt)
%DATA_CUBE: Computes a seismic cube composed of the superpostion of multidimensional plane waves.
%           This code is used to test reconstruction algorithms. 
%
% [d] = data_cube(N,dt,f0,nt,dx,t0,p,A,opt);
% 
% IN    N:      spatial dimensions of cube containing ne events (Plane waves) N(1:nsd)
%               where nsd=number of spatial dimensions 
%       dt:     sampling interval in seconds 
%       f0:     central freq. of wavelet
%       nt:     number of time samples
%       dx:     sampling increment in each spatial dimension of the cube  
%       t0:     intercept  t0(1:ne) 
%       p:      ray-parameter  p(1:nsd, 1:ne) 
%       A:      Amplitudes A(1:ne) 
%       option: 'linear' or 'parabolic'
%
% 
% OUT   d:      seismic cube  d(1:nt, 1:n1) :                        2D when N=[n1]
%               seismic cube  d(1:nt, 1:n1, 1:n2) :                  3D when N=[n1,n2]
%               seismic cube  d(1:nt, 1:n1, 1:n2, 1:n3) :            4D when N=[n1,n2,n3]
%               seismic cube  d(1:nt, 1:n1, 1:n2, 1:n3, 1:n4) :      5D when N=[n1,n2,n3,n4]
%
% 
% Notes:
% 
%       For option 'parabolic', p is residual moveout at far offet 
% 
%       data_cube() produce a 3D cube with default parameters. 
% 
% Examples: 
%
%  1) data in t-x (2D cube) with 2 linear events 
% 
%     N = 10;
%     nt = 200;     
%     dx = 10.; 
%     dt = 4./1000.;
%     p  = [0.2, -0.1]
%     t0 = [ 0.1,  0.4];
%     A  = [-1.0,  1.0];
%     f0 = 10.;
%     opt = 'parabolic';
%     [d] = data_cube(N,dt,f0,nt,dx,t0,p,A,opt);
%     wigb(d); 
%
%  2) data in t-x1-x2 (3D cube) with 2 linear events 
%     N = [40,20] 
%     nt = 140;     
%     dx = [10.,10.]; 
%     dt = 4./1000.;
%     p  = [0.0001, -0.0001;
%           0.0002, -0.0001]
%     t0 = [0.1, 0.2];
%     A  = [-1.0, 1.0];
%     f0 = 10.;
%     opt = 'linear'
%     [d] = data_cube(N,dt,f0,nt,dx,t0,p,A,opt);
%     imagesc(reshape(d,nt,N(1)*N(2)));
%
%  3) data in t-x1-x2-x3-x4  (5D cube) with 2 linear events 
%     N = [10,14,13,50] 
%     nt = 200;     
%     dx = [10.,10.,10.,10.]; 
%     dt = 4./1000.;
%     p  = [0.0200, -0.1000;
%           0.0300,  0.0000;
%           0.0100,  0.0100;
%           0.0020, -0.0200]
%     t0 = [0.1, 0.2];
%     A  = [-1.0, 1.0];
%     f0 = 10.;
%     opt = 'parabolic'
%     [d] = data_cube(N,dt,f0,nt,dx,t0,p,A,opt);
%     imagesc(reshape(d,nt,N(1)*N(2)*N(3)*N(4)));
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


if nargin == 0
     N = [40,20] 
     nt = 140;     
     dx = [10.,10.]; 
     dt = 4./1000.;
     p  = [0.0001, -0.0001;
           0.0002, -0.0001]
     t0 = [0.1, 0.2];
     A  = [-1.0, 1.0];
     f0 = 10.;
     opt = 'linear'
     [d] = data_cube(N,dt,f0,nt,dx,t0,p,A,opt);
end;



ND = length(N);


%--------------- 2D case --------------------

if ND == 1 

 n1 = N(1);
 ne = length(p);
 nfft = 2*(2^nextpow2(nt));

 x1 = [0:1:n1-1]*dx(1);

if strcmp(opt,'parabolic')
  x1 = (x1/max(x1)).^2;  
end

 i = sqrt(-1);
 D = zeros(nfft,n1);
 [wavelet,tw] = ricker(f0,dt);
 nwavelet = length(wavelet);
 tdelay = nwavelet*dt/2;
 W = fft(wavelet,nfft);
  for iw=2:nfft/2;
     w = 2*pi*(iw-1)/(nfft*dt);
     T = zeros(1,n1);
      for ie=1:ne;
        x11 = exp(-i*w*x1*p(1,ie));
        T = T + A(ie)*exp(-i*w*(t0(ie)-tdelay))*W(iw)*x11;
     end;
   D(iw,:) =  T;
   D(nfft-iw+2,:) = conj(T);
  end;

 d = ifft(D,[],1);
 d = real(d(1:nt,:));

end;

%--------------- 3D case --------------------

if ND == 2 

 n1 = N(1);
 n2 = N(2);
 ne = size(p,2);
 nfft = 2*(2^nextpow2(nt));

 x1 = [0:1:n1-1]*dx(1);
 x2 = [0:1:n2-1]*dx(2);


if strcmp(opt,'parabolic')
  x1 = (x1/max(x1)).^2;  
  x2 = (x2/max(x2)).^2;  
end

 i = sqrt(-1);
 D = zeros(nfft,n1,n2);
 [wavelet,tw] = ricker(f0,dt);
 nwavelet = length(wavelet);
 tdelay = nwavelet*dt/2;
 W = fft(wavelet,nfft);
  for iw=2:nfft/2;
     w = 2*pi*(iw-1)/(nfft*dt);
     T = zeros(n1,n2);
      for ie=1:ne;
        x11 = exp(-i*w*x1*p(1,ie));
        x22 = exp(-i*w*x2*p(2,ie));
        [M1,M2] = ndgrid(x11,x22);
        T = T + A(ie)*exp(-i*w*(t0(ie)-tdelay))*W(iw)*M1.*M2;
     end;
   D(iw,:,:) =  T;
   D(nfft-iw+2,:,:) = conj(T);
  end;

 d = ifft(D,[],1);
 d = real(d(1:nt,:,:));

end;


%--------------- 4D case --------------------

if ND == 3 

 n1 = N(1);
 n2 = N(2);
 n3 = N(3);
 ne = size(p,2);
 nfft = 2*(2^nextpow2(nt));

 x1 = [0:1:n1-1]*dx(1);
 x2 = [0:1:n2-1]*dx(2);
 x3 = [0:1:n3-1]*dx(3);

if strcmp(opt,'parabolic')
  x1 = (x1/max(x1)).^2;  
  x2 = (x2/max(x2)).^2;  
  x3 = (x3/max(x3)).^2;  
end

 i = sqrt(-1);
 D = zeros(nfft,n1,n2,n3);
 [wavelet,tw] = ricker(f0,dt);
 nwavelet = length(wavelet);
 tdelay = nwavelet*dt/2;
 W = fft(wavelet,nfft);
  for iw=2:nfft/2;
     w = 2*pi*(iw-1)/(nfft*dt);
     T = zeros(n1,n2,n3);
      for ie=1:ne;
        x11 = exp(-i*w*x1*p(1,ie));
        x22 = exp(-i*w*x2*p(2,ie));
        x33 = exp(-i*w*x3*p(3,ie));
        [M1,M2,M3] = ndgrid(x11,x22,x33);
        T = T + A(ie)*exp(-i*w*(t0(ie)-tdelay))*W(iw)*M1.*M2.*M3;
     end;
   D(iw,:,:,:) =  T;
   D(nfft-iw+2,:,:,:) = conj(T);
  end;

 d = ifft(D,[],1);
 d = real(d(1:nt,:,:,:));

end;

%--------------- 5D case --------------------

if ND == 4 

 n1 = N(1);
 n2 = N(2);
 n3 = N(3);
 n4 = N(4);
 ne = size(p,2);
 nfft = 2*(2^nextpow2(nt));

 x1 = [0:1:n1-1]*dx(1);
 x2 = [0:1:n2-1]*dx(2);
 x3 = [0:1:n3-1]*dx(3);
 x4 = [0:1:n4-1]*dx(4);

if strcmp(opt,'parabolic')
  x1 = (x1/max(x1)).^2;  
  x2 = (x2/max(x2)).^2;  
  x3 = (x3/max(x3)).^2;  
  x4 = (x4/max(x4)).^2;  
end

 i = sqrt(-1);
 D = zeros(nfft,n1,n2,n3,n4);
 [wavelet,tw] = ricker(f0,dt);
 nwavelet = length(wavelet);
 tdelay = nwavelet*dt/2;
 W = fft(wavelet,nfft);
  for iw=2:nfft/2;
     w = 2*pi*(iw-1)/(nfft*dt);
     T = zeros(n1,n2,n3,n4);
      for ie=1:ne;
        x11 = exp(-i*w*x1*p(1,ie));
        x22 = exp(-i*w*x2*p(2,ie));
        x33 = exp(-i*w*x3*p(3,ie));
        x44 = exp(-i*w*x4*p(4,ie));
        [M1,M2,M3,M4] = ndgrid(x11,x22,x33,x44);
        T = T + A(ie)*exp(-i*w*(t0(ie)-tdelay))*W(iw)*M1.*M2.*M3.*M4;
     end;
   D(iw,:,:,:,:) =  T;
   D(nfft-iw+2,:,:,:,:) = conj(T);
  end;

 d = ifft(D,[],1);
 d = real(d(1:nt,:,:,:,:));

end

return;
