function [df] = fx_decon(d,dt,lf,mu,flow,fhigh);
%FX_DECON: SNR enhancement using fx-deconvolution.
%
%  [DATA_f] = fx_decon(D1,D2,dt,lf,mu,flow,fhigh);
% 
%  IN   d:      the data (a matrix, traces are columns) 
%       dt:     sampling interval in sec
%       lf:     length of operator (length of the filter)
%       mu:     pre-whitening 
%       flow:   min  freq. in the data in Hz
%       fhigh:  max  freq. in the data in Hz
% 
%  OUT  df:     filtered data 
%
% 
%  Reference: Canales, 1984, Random noise reduction, 54.th. Ann. Internat. 
%             Mtg., Soc. Expl. Geophys., Expanded Abstracts, pp. 525-527
%
%  Note: Canales' method is modified to use non Toeplitz system of equations
%        with backward and forward prediction filters
%       
%  Example: see fx_decon_demo.m
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


 [nt,ntraces] = size(d);
 nf = 2^nextpow2(nt);
 
 D_f = zeros(nf,ntraces);
 D_b = zeros(nf,ntraces);

% First and last samples of the DFT.

 ilow  = floor(flow*dt*nf)+1; 

  if ilow<1; 
   ilow=1; 
  end;

 ihigh = floor(fhigh*dt*nf)+1;

  if ihigh > floor(nf/2)+1; 
   ihigh=floor(nf/2)+1; 
  end

% Transform to FX

 D = fft(d,nf,1);

 for k = ilow:ihigh;
  aux  = D(k,:)';
  [aux_out_f,aux_out_b] = ar_modeling(aux,lf,mu);
  D_f(k,:) = aux_out_f';
  D_b(k,:) = aux_out_b';
 end;

% Honor symmetries

 for k=nf/2+2:nf
  D_f(k,:) = conj(D_f(nf-k+2,:));
  D_b(k,:) = conj(D_b(nf-k+2,:));
 end

% Back to TX (the output) 

 d_f = real(ifft(D_f,[],1));
 d_f = d_f(1:nt,:);

 d_b = real(ifft(D_b,[],1));
 d_b = d_b(1:nt,:);

% Average predictions (forward and backward)

 df = d_f+d_b;
 df(:,lf+1:ntraces-lf)= df(:,lf+1:ntraces-lf)/2;
 
return

function [yf,yb] = ar_modeling(x,lf,mu);
%AR_MODELING: autoregressive modeling of 1D spatial data
%
%  IN    x:   data 
%        lf:  length of the operator
%        mu:  pre-whitening in %
%      
%  OUT   yf:  prediction of the data using forward AR modelling
%        yb:  prediction of the data using backward AR modelling
% 
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.



   nx = length(x);

% backward AR modelling

   y  = x(1:nx-lf,1);
   C  = x(2:nx-lf+1,1);
   R  = x(nx-lf+1:nx,1);
   M = hankel(C,R);

   B = M'*M;  
   beta = B(1,1)*mu/100;
   ab = (B + beta*eye(lf))\M'*y;

   temp = M*ab;
   temp = [temp;zeros(lf,1)];
   yb = temp;


% forward  AR modelling

   y  = x(lf+1:nx,1);
   C  = x(lf:nx-1,1);
   R = flipud(x(1:lf,1));
   M = toeplitz(C,R);
   B = M'*M;  
   beta = B(1,1)*mu/100;
   af = (B + beta*eye(lf))\M'*y;
   temp = M*af;
   temp = [zeros(lf,1);temp];
   yf = temp;

return
