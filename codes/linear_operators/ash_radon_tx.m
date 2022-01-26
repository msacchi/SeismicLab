function [Out] = ash_radon_tx(In, Par, itype);
%RADON_TX: Operators for t-x, tau-v-a forward and adjoint Radon transforms
%          to compute  the apex shifted hyperbolic Radon transform. 
%
%  [Out] = asr_radon_tx(In, Par);
%
%  IN    Radon coefficients if itype =  1 (Forward transform)  with In(nt,nv,na)
%        CSG or CRG gather  if itype = -1 (Adjoint transform)  with In(nt,nh)
%
%  OUT   CSG or CRG gather  if itype =  1 (Forward transform)  with Out(nt,nh)
%        Radon coeffcients  if itype = -1 (Adjoint transform)  with Out(nt,nv,na) 
%  
%  Par.h       :  vector containing the nh offsets 
%  Par.v       :  vector containing the nv velocities in m/s 
%  Par.a       :  vector containing the na apexes in m
%  Par.dt      :  sampling interval 
%  Par.f       :  frequency corners of BP operator - this acts like a zero phase wavelet. 
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
   v = Par.v;
   a = Par.a;
  dt = Par.dt;
  f  = Par.f;  

  f1 = f(1);
  f2 = f(2);
  f3 = f(3);
  f4 = f(4);

 nv = length(v);
 na = length(a);
 nh = length(h);

 if itype == 1; 
    m = In;
    [nt, nv, na] = size(m);
    d = zeros(nt,nh);
    m = bp_filter(m,dt,f1,f2,f3,f4);
 end;

 if itype == -1; 
    d = In;
    [nt, nh] = size(d);
    m = zeros(nt,nv,na);
 end;

 time  = @(x1,x2,x3,x4) sqrt(x1^2 - (x2-x3)^2/x4^2); 

    for it = 1:nt 
     for iv = 1:nv
      for ia = 1:na
       for ih = 1:nh 

       tau = time((it-1)*dt, h(ih), a(ia), v(iv));

       if tau >= 0.;  

         itau  = tau/dt+1;
         itau1 = floor(itau);
         itau2 = itau1 + 1;
         alpha = itau-itau1;

         if itau2<nt & itau1>1;

          if itype == 1; 
                 d(it,ih) = d(it,ih) + (1-alpha) * m(itau1,iv,ia) + alpha*m(itau2,iv,ia);
          else 
                 m(itau1,iv,ia) = m(itau1,iv,ia) + (1-alpha) * d(it,ih);
                 m(itau2,iv,ia) = m(itau2,iv,ia) +     alpha * d(it,ih);
         end
      end;

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
