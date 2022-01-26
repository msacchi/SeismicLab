function [dout,e1,e2,freq] = pocs(d,dtrue,T,dt,f_low,f_high, option, perc_i, perc_f, N, a, tol);
%POCS: 3D (x,y,t) seismic data regularization using POCS
%
% [dout,converg] = pocs(d,dtrue,T,dt,f_low,f_high, option, perc_i, perc_f, N, a, tol);
%
%  IN   d:      data  a t,x,y cube)
%       dtrue:  true data prior to decimation (insert d is you do not have dtrue)
%       T:      sampling operator
%       f_low:  min freq in Hz in the data
%       f_high: max freq in Hz inthe data
%       option:    = 1 linear threshold schedulde
%                    2 exponential threshold schedulde
%                    3 data adaptive threshold schedule
%       perc_i: percentage of max amplitude for initial threshold   
%       perc_f: percentage of max amplitude for final threshold (perc_f<<perc_i)   
%       N:      maximum number of iterations 
%       a:      a = 1 for data with high SNR or when one wants to fit the noise 
%               a < 1 for noise attenuation (try a=0.2 to 0.4)
%       tol:    tolerance (try 0.01). This is to truncate the number of iterations 
%
%  OUT  dout:  reconstructed cube
%       e1:    relative decrease in cost versus iteration and freq
%       e2:    relative decrease in mse versus iteration and freq (can only be computed if
%              dtrue is provided)
%              e1 and e2 are used to analyze convergence 
%

 [nt,nx1,nx2] = size(d);

  nf = 2^nextpow2(nt);
 nk1 = 2^nextpow2(nx1);
 nk2 = 2^nextpow2(nx2);

% T is sample, S is the reinsertion operator 
 
 S = 1.-T;

% Min and max freq indeces for the reconstruction

 k_low = max(floor(f_low*dt*nf) + 1,2);
 k_high = min(floor(f_high*dt*nf) + 1,nf/2);
 
 Dout = zeros(nf,nx1,nx2);
 D = fft(d,nf,1);
 Dtrue = fft(dtrue,nf,1);
 
 for k = k_low:k_high;

  freq(k) = (k-1)/(nf*dt);
  x = squeeze(D(k,:,:));
  xtrue = squeeze(Dtrue(k,:,:));
  th = th_schedule(option,x,perc_i,perc_f,N);
  y = x;

   iter = 1;
   E1 = 100;

  while (iter <=N & E1>tol);

   Y = fft2(y,nk1,nk2);
   A = abs(Y); Angle = angle(Y);
   I = find(A<th(iter)); A(I)=0; Y = A.*exp(i*Angle);
   y = ifft2(Y); y = y(1:nx1,1:nx2);
   yold = y;
   y = a*x + (1-a)*T.*y + S.*y;
   dif1 = y-yold;
   dif2 = y-xtrue;
   c1 = sum( (abs(dif1(:))).^2);
   c2 = sum( (abs(dif2(:))).^2);
   c  = sum( (abs(y(:))).^2);
   ct = sum( (abs(xtrue(:))).^2);
   E1 = c1/c;
   E2 = c2/ct;
   e1(k-k_low+1,iter)=E1;
   e2(k-k_low+1,iter)=E2;
   iter = iter + 1;

  end

  iter-1;
   
  Dout(k,:,:) = y;
  Dout(nf-k+2,:,:)  = conj(y);

 end;

 dout = ifft(Dout,[],1);
 dout = dout(1:nt,:,:);

 return

function th = th_schedule(option,x,perc_i,perc_f,N);
% Function used by pocs to define the threshold schedule based 
% on parameter option
%
%  option == 1 --> linear
%  option == 2 --> exponential
%  option == 3 --> use real amplitude to define curve
%

  [nx1,nx2] = size(x);
  nk1 = 2^nextpow2(nx1);
  nk2 = 2^nextpow2(nx2);
  X = fft2(x,nk1,nk2);
  A = abs(X); 
  Amax = max(A(:));
  th_i = perc_i*Amax/100;
  th_f = perc_f*Amax/100;
  k = [1:1:N];

 if option==1;
  th = th_i + (th_f-th_i)*(k-1)/(N-1);
 end;
  
 if option==2;
  b = -log(th_f/th_i);
  th = th_i*exp(-b*(k-1)/(N-1));
 end;

 if option==3 ; 
  Amax = max(A(:));
  avec=sort(reshape(A,nk1*nk2,1),'descend');
  I = find(avec>perc_i*Amax/100);
  avec(I)=[];
  I = find(avec<perc_f*Amax/100);
  avec(I)=[];
  th = zeros(N,1);
  th(1) = avec(1);
  th(N) = avec(end);
  for j=2:N-1
  th(j)=avec(ceil((j-1)*length(avec)/(N-1)));
  end
 end;

 return;
