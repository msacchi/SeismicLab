function [dout] = mwni(data_binned, T, dt, fmin, fmax, NDim, K, max_iter_int, max_iter_ext)
% MWNI: reconstruction of binned data using Minimum Weigthed Norm Interpolation
% 
% References: 
%
%      Minimum weighted norm interpolation of seismic records
%      Bin Liu, Mauricio D. Sacchi
%      GEOPHYSICS Nov 2004, Vol. 69, No. 6, pp. 1560-1568
%
%      Five-dimensional interpolation: Recovering from acquisition constraints
%      Daniel Trad
%      GEOPHYSICS Nov 2009, Vol. 74, No. 6, pp. V123-V132
%
%
% IN.  data_binned is the data binned with zeros in missing positions 
%      data_binned(nt,n1)                   if NDim = 2
%      data_binned(nt,n1,n2)                if NDim = 3
%      data_binned(nt,n1,n2,n3)             if NDim = 4
%      data_binned(nt,n1,n2,n3,n4)          if NDim = 5
%  
%      T is the sampling operator with 1 when the trace is alive and 0 when the trace is dead 
%      T(n1)                                if NDim = 2
%      T(n1,n2)                             if NDim = 3
%      T(n1,n2,n3)                          if NDim = 4
%      T(n1,n2,n3,n4)                       if NDim = 5
%      
%      dt:   is sampling interval in secs
%      fmin: is minimum freq in Hz to process
%      fmax: is maximum freq in Hz to process
%      NDim: number of dimensions of data_binned including time (2,3,4,5)
%      K:    number of wavenumbers in each direction 
%      K = k1              if NDim = 2  k1>n1
%      K = [k1,k2]         if NDim = 3  k1>n2 and k2>n2
%      K = [k1,k2,k3]      if NDim = 4
%      K = [k1,k2,k3,k4]   if NDim = 5
%      iter_max_int: max number of internal iterations for IRLS
%      iter_max_ext: max number of externa iterations for IRLS
%
% OUT. dout is the reconstruced data of the same size of data_binned but with dead traces reconstructed
%
% ------------------------------------------
% Mauricio D Sacchi, Natal, August 2015
% msacchi@ualberta.ca
% ------------------------------------------
%

  if NDim == 1; error('Your data needs to have at least 2 dimensions, time and one spatial coordinate');end
  if NDim == 2;  [nt,n1]              = size(data_binned); N = n1; end
  if NDim == 3;  [nt,n1,n2]           = size(data_binned); N = [n1,n2]; end
  if NDim == 4;  [nt,n1,n2,n3]        = size(data_binned); N = [n1,n2,n3]; end
  if NDim == 5;  [nt,n1,n2,n3,n4]     = size(data_binned); N = [n1,n2,n3,n4];end

% Number of temporal frequencies of the DFT (20% padding)

  nfft = 2*floor(0.5*nt*1.2)              

% Evaluate 1D DFT along time for all traces in the volume data_binned

  D = fft(data_binned,nfft,1); 

% Indeces of max and min frequencies for the DFT

  ifreq_min = floor(fmin*dt*nfft) + 1;
  ifreq_max = floor(fmax*dt*nfft) + 1;

% Set initial volumes according to NDim

  if NDim == 2;  k1 = K(1);                                 x0 = zeros(k1,1);           dout = zeros(nfft,n1); end;
  if NDim == 3;  k1 = K(1); k2 = K(2);                      x0 = zeros(k1,k2);          dout = zeros(nfft,n1,n2); end;
  if NDim == 4;  k1 = K(1); k2 = K(2); k3 = K(3);           x0 = zeros(k1,k2,k3);       dout = zeros(nfft,n1,n2,n3); end;
  if NDim == 5;  k1 = K(1); k2 = K(2); k3 = K(3); k4=K(4);  x0 = zeros(k1,k2,k3,k4);    dout = zeros(nfft,n1,n2,n3,n4); end;


% Solve problem in space one frequency at a time using IRLS (see paper by Liu and Sacchi, 2004 to understand IRLS)
% Aaccording to NDim I can handle different type of interpolation 2D,3D,4D,5D... notice that in all cases
% we first estimate the sparse Fourier wavenumber coefficients with mwni_irls and then we use operator_nfft
% with unity sampling to recover the data.

% ------------------------------------------
% ------------    2D case  -----------------
% ------------------------------------------


 if NDim == 2; 
  for ifreq = ifreq_min:ifreq_max;
    y = squeeze(D(k,:));
    [x,misfit] = mwni_irls(y,x0,T,NDim,M,K,itmax_internal,itmax_external,silence);
    y = operator_nfft(x,ones(n1,1),NDim,N,K,1);
    Dout(k,:) = y;
  end
 end;

% ------------------------------------------
% ------------    3D case  -----------------
% ------------------------------------------

 if NDim == 3; 
  for ifreq = ifreq_min:ifreq_max;
    y = squeeze(D(k,:,:));
    [x,misfit] = mwni_irls(y,x0,T,NDim,M,K,itmax_internal,itmax_external,silence);
    y = operator_nfft(x,ones(n1,n2),NDim,N,K,1);
    Dout(k,:,:) = y;
  end
 end;

% ------------------------------------------
% ------------    4D case  -----------------
% ------------------------------------------

 if NDim == 4; 
  for ifreq = ifreq_min:ifreq_max;
    y = squeeze(D(k,:,:,:));
    [x,misfit] = mwni_irls(y,x0,T,NDim,M,K,max_iter_int,max_iter_ext,silence);
    y = operator_nfft(x,ones(n1,n2,n3),NDim,N,K,1);
    Dout(k,:,:,:) = y;
  end
 end;

% ------------------------------------------
% ------------    5D case  -----------------
% ------------------------------------------

 if NDim == 5; 
  for ifreq = ifreq_min:ifreq_max;
    ifreq
    y = squeeze(D(ifreq,:,:,:,:));
    [x,misfit] = mwni_irls(y,x0,T,NDim,N,K,max_iter_int,max_iter_ext,1);
    y = operator_nfft(x,ones(n1,n2,n3,n4),NDim,N,K,1);
    dout(ifreq,:,:,:,:) = y;
  end
 end;

% Fourier domain symmetries  to guarantee that time series are real when we return to time 

 if NDim == 2; 
  for k = nfft/2+2:nfft;
   Dout(k,:) = conj(Dout(nfft-k+2,:));
  end;
  dout = ifft(Dout,nfft,1); 
  dout = real(dout(1:nt,:)); 
 end

 if NDim == 3; 
  for k = nfft/2+2:nfft;
   Dout(k,:,:) = conj(Dout(nfft-k+2,:,:));
  end;
  dout = ifft(Dout,nfft,1); 
  dout = real(dout(1:nt,:,:)); 
 end

 if NDim == 4; 
  for k = nfft/2+2:nfft;
   Dout(k,:,:,:) = conj(Dout(nfft-k+2,:,:,:));
  end;
  dout = ifft(Dout,nfft,1); 
  dout = real(dout(1:nt,:,:,:)); 
 end


 if NDim == 5; 
  for k = nfft/2+2:nfft;
   dout(k,:,:,:,:) = conj(dout(nfft-k+2,:,:,:,:));
  end;
  dout = ifft(dout,nfft,1);   % 1D ifft  
  dout = real(dout(1:nt,:,:,:,:)); 
 end



return
