function [Out] = radon_fx(In, Par, itype);
%RADON_FX: Operators for f-x, tau-p forward and adjoint Radon transforms
%          to compute linear or  parabolic Radon transform. 
%
%  [Out] = radon_forward_adjoint(In, Par, type_of_transform);
%
%  IN    Radon coefficients if itype =  1 (Forward transform)  with In(np,nt)
%        CMP gather         if itype = -1 (Adjoint transform)  with In(nh,nt)
%
%  OUT   CMP gather         if itype =  1 (Forward transform)  with Out(nh,nt)
%        Radon coeffcients  if itype = -1 (Adjoint transform)  with Out(np,nt) 
%  
%  Par.h       :  vector containing the nh offsets 
%  Par.p       :  vector containing the np curvatures normalized by far offset (s) if itype='parab'
%  Par.p       :  vector containing the np dips in s/m if itype='linear'
%  Par.dt      :  sampling interval 
%  Par.f       :  frequency corners of BP operator - this acts like a zero phase wavelet. 
%  Par.transf  :  'parab', 'linear', 'hyperb'
%
%  itype =  1  :  forward
%  itype = -1  :  adjoint
%
%  Notes: 
% 
%  This function calls bp_filter.m from SeismicLab. It is like doing the deconvoluted RT with
%  a band-pass zero phase wavelet. 
%
%
% M D Sacchi
% msacchi@ualberta.ca 
%


   dt = Par.dt;
    h = Par.h;
    p = Par.p;
    f = Par.f;
    
   f1 = f(1);
   f2 = f(2);
   f3 = f(3);
   f4 = f(4);
  
if  strcmp(Par.transf,'parab') ; hmax = max(abs(h)); h=(h/hmax).^2;   end; 

nh = length(h);
np = length(p);

if itype == 1;
    nt = size(In,1); 
    m = In;
    m = bp_filter(m,dt,f1,f2,f3,f4);
else
    nt = size(In,1); 
    d = In;
end

 nfft = 2*2^nextpow2(nt);
 ilow  = floor(f(1)*dt*nfft)+1;
 ihigh = floor(f(4)*dt*nfft)+1;

 i = sqrt(-1);
 A = h'*p;

if itype == 1;
    
    M = fft(m,nfft,1);
    D = zeros(nfft,nh);
    
   for ifreq = ilow:ihigh
        w = 2.*pi*(ifreq-1)/nfft/dt;
        L = exp(i*w*A);
        x = M(ifreq,:)';
        y = L * x;
        D(ifreq,:) = y';
    end
    
    for k=nfft/2+2:nfft
        D(k,:) = conj(D(nfft-k+2,:));
    end
    
    d = ifft(D,[],1);
    d = real( d(1:nt,:));
    
    Out = d;
    
else
    
    D = fft(d,nfft,1);
    M = zeros(nfft,np);
    
    for ifreq = ilow:ihigh
        w = 2.*pi*(ifreq-1)/nfft/dt;
        L = exp(i*w*A);
        y = D(ifreq,:)';
        x = L' * y;
        M(ifreq,:) = x';
    end
    
    for k = nfft/2+2:nfft
        M(k,:) = conj(M(nfft-k+2,:));
    end
    
    m = ifft(M,[],1);
    m = real( m(1:nt,:));
    
    Out = m;
    
end

if itype == -1;
    m = bp_filter(m,dt,f1,f2,f3,f4);
    Out = m;
end;


return;



