function [Out] = radon_tx(In, Par, itype);
%RADON_TX: Operators for t-x, tau-v forward and adjoint Radon transforms
%          to compute linear, parabolic or hyperbolic Radon transform. 
%
%  [Out] = radon_forward_adjoint(In, Par, type_of_transform);
%
%  IN    Radon coefficients if itype =  1 (Forward transform)  with In(nv,nt)
%        CMP gather         if itype = -1 (Adjoint transform)  with In(nh,nt)
%
%  OUT   CMP gather         if itype =  1 (Forward transform)  with Out(nh,nt)
%        Radon coeffcients  if itype = -1 (Adjoint transform)  with Out(nv,nt) 
%  
%  Par.h       :  vector containing the nh offsets 
%  Par.p       :  vector containing the np velocities in m/s  if itype='hyperb' 
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

   h = Par.h;
   p = Par.p;
  dt = Par.dt;
  f  = Par.f;  

  f1 = f(1);
  f2 = f(2);
  f3 = f(3);
  f4 = f(4);

 np = length(p);
 nh = length(h);

 if itype == 1; 
    m = In;
    [nt, np] = size(m);
    d = zeros(nt,nh);
    m = bp_filter(m,dt,f1,f2,f3,f4);
 end;

 if itype == -1; 
    d = In;
    [nt, nh] = size(d);
    m = zeros(nt,np);
 end;

if  strcmp(Par.transf,'parab') ; hmax = max(abs(h));  time  = @(a,b,c) a + (b/hmax)^2*c; end; 
if  strcmp(Par.transf,'linear');                      time  = @(a,b,c) a + (b*c); end; 
if  strcmp(Par.transf,'hyperb');                      time  = @(a,b,c) sqrt(a^2 + (b^2/c^2)); end; 

 tmax = (nt-1)*dt; 
 tau = [0:1:nt-1]*dt; 

if itype == 1; 

  for itau = 1:1:nt 
    for ih = 1:nh 
     for ip = 1:np 
      t = time(tau(itau), h(ih), p(ip));
       if t>0 & t<tmax 
         it    = t/dt+1;
         it1   = floor(it);
         it2   = it1 + 1;
         alpha = it  - it1;
                 d(it1,ih) = d(it1,ih) + (1-alpha) * m(itau,ip);
                 d(it2,ih) = d(it2,ih) +     alpha * m(itau,ip);
       end
     end
    end
   end
 
 else 

  for itau = 1:1:nt 
    for ih = 1:nh 
     for ip = 1:np 
      t = time(tau(itau), h(ih), p(ip));
       if t>0 & t<tmax 
         it    = t/dt+1;
         it1   = floor(it);
         it2   = it1 + 1;
         alpha = it  - it1;
                 m(itau,ip) = m(itau,ip) + (1-alpha) * d(it1,ih) + alpha*d(it2,ih);
       end
     end
    end
   end

 end 


 if itype == 1; 
    Out = d;
 end;

 if itype == -1; 
    m = bp_filter(m,dt,f1,f2,f3,f4);
    Out = m;
 end;

 return
