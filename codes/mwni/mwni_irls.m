function [x,misfit] = mwni_irls(d,m0,T,NDim,N,K,itmax_internal,itmax_external,silence);
% MWNI_IRLS: minimize   ||A x - d||_2 with x sparse where x are space domain Fourier coefficients
%            and d is spatial data for one monocromatic frequency.  
%            IRLS -> Iterative Reweigthed Least Squares 
% 
% References: 
%
%      Minimum weighted norm interpolation of seismic records
%      Bin Liu, Mauricio D. Sacchi
%      GEOPHYSICS Nov 2004, Vol. 69, No. 6, pp. 1560-1568
%
% IN.  d complex spatial data
%      d(n1)                   if NDim = 2
%      dn1,n2)                 if NDim = 3
%      d(n1,n2,n3)             if NDim = 4
%      d(n1,n2,n3,n4)          if NDim = 5
%  
%      T is the sampling operator with 1 when the trace is alive and 0 when the trace is dead 
%      T(n1)                                if NDim = 2
%      T(n1,n2)                             if NDim = 3
%      T(n1,n2,n3)                          if NDim = 4
%      T(n1,n2,n3,n4)                       if NDim = 5
%      
%      NDim: number of dimensions of data_binned including time (2,3,4,5)
%      K:    number of wavenumbers in each direction 
%      K = k1              if NDim = 2
%      K = [k1,k2]         if NDim = 3
%      K = [k1,k2,k3]      if NDim = 4
%      K = [k1,k2,k3,k4]   if NDim = 5
%      iter_max_int: max number of internal iterations for IRLS
%      iter_max_ext: max number of externa iterations for IRLS
%      silence = 1 (silent mode), 0 (will report Misfit versus Iteration)
%
% OUT. x are the inverted spatial wavenumber Fourier coefficients
%
% ------------------------------------------
% Mauricio D Sacchi, Natal, August 2015
% msacchi@ualberta.ca
% ------------------------------------------


  z = zeros(size(m0));
  P = ones(size(z));
  kc = 1;
  Misfit = [];
  x = z;
  for l = 1:itmax_external
  di = operator_nfft(P.*z,T,NDim,N,K,1);

  r = d-di;
  g = operator_nfft(r,T,NDim,N,K,-1);
  g = g.*P;
  s = g;
  gammam = cgdot(g);
  k = 1;
  while  k<=itmax_internal;

   ss = operator_nfft(P.*s,T,NDim,N,K,1);
   den = cgdot(ss);
   alpha = gammam/(den+1.e-8);
   z = z+alpha*s;
   r = r-alpha*ss;
   misfit(kc) =  cgdot(r);
   g = operator_nfft(r,T,NDim,N,K,-1);
   g = g.*P;
   gamma = cgdot(g);
   beta = gamma/(gammam + 1.e-7);
   gammam = gamma;
   s = g + beta*s;
 if silence == 0;
   say = sprintf(['Iteration = %d Misfit=%0.5g'],k,misfit(k));
   disp(say);
 end;
   k = k + 1;
   kc = kc + 1;
  end;

  x = P.*z;
  y = x/max(abs(x(:)));
  P = abs(y)+0.001;
   
 end;
 return


function [out] = cgdot(in);

% Dot product
   
 temp  =   in.*conj(in);
 out = sum(temp(:));

return;

