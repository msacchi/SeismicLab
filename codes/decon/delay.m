function dout = delay(d1,d2,max_delay);
%DELAY: Delay d2 toward d1
%
%  [dout] = delay(d1,d2,max_delay);
%
%  IN   d1:        Reference data (vector or matrix, traces are columns)
%       d2:        Data to delay (d1 and d2 must have the same size)
%       max_delay: Max x-correlation lag
%
%  OUT  dout:      Delayed version of d2 toward reference data d1
%
%  Note: Use this program to apply a +/- time shift after deconvolution.
%
%        
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.
%

  [nt,nx] = size(d1);

  r = zeros(2*max_delay+1,1);

   for k = 1:nx
    temp = xcorr(d1(:,k),d2(:,k),max_delay);
    r = r + temp;
   end;

 [rmax,Lag] = max(r);
 time_to_move = Lag-max_delay-1

 if time_to_move>0;
  dout = [zeros(time_to_move,nx);d2];
  dout = dout(1:nt,:);
 end;
    
 if time_to_move ==0; dout=d2; end;

 if time_to_move<0; 
  ll = -time_to_move+1;
  dout = [d2(ll:nt,:);zeros(ll-1,nx)];
 end;
 
 return;
