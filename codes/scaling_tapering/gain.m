function [dout] = gain(d,dt,option1,parameters,option2);
%GAIN: Gain a group of traces.
%
%  [dout] = gain(d,dt,option1,parameters,option2);
%
%  IN   d:         traces, a matrix, a cube etc, first dimension is time 
%       dt:        ampling interval
%       option1 = 'time' parameters = [a,b],  gain = t.^a . * exp(-bt)
%               = 'agc' parameters = [agc_gate], length of the agc gate in secs
%       option2 = 0  No normalization
%               = 1  Normalize each trace by amplitude
%               = 2  Normalize each trace by rms value
%
%  OUT  dout(nt,nx): traces after application of gain function
%
%
%  Example:
%
%    d = hyperbolic_events; dout = gain(d,0.004,'agc',0.05,1);
%    wigb([d,dout]);
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


 N     = size(d);
 nt    = N(1);
 x     = reshape(d,nt,prod(N(2:end)));

 if strcmp(option1,'time');   % Geometrical spreading-like gain

  a = parameters(1);
  b = parameters(2);
  t = (0:1:nt-1)*dt; 
  tgain = (t.^a).*exp(b.*t);

  for k = 1: prod(N(2:end));
   xout(:,k)  = x(:,k).*tgain';
  end;

 end

 if strcmp(option1,'agc');    % AGC 

  L = parameters(1)/dt+1;
  L = floor(L/2);
  h = hamming(2*L+1);

  for k = 1:prod(N(2:end));
   aux =  x(:,k);
   e = aux.^2;
   rms = sqrt(conv2(e,h,'same'));
   epsi = 1.e-10*max(rms);
   op = rms./(rms.^2+epsi);
   xout(:,k) = x(:,k).*op;
   end
  end

 if option2==1;                % Normalize by amplitude 

   for k = 1:prod(N(2:end));
    aux =  xout(:,k);
    amax = max(abs(aux));
    xout(:,k) = xout(:,k)/(amax+0.00000001);
   end

 end


 if option2==2;                % Normalize by rms 

   for k = 1:prod(N(2:end));
    aux =  xout(:,k);
    amax =  sqrt(sum(aux.^2)/nt);
    xout(:,k) = xout(:,k)/(amax+0.00000001);
   end

 end

dout = reshape(xout,[nt,N(2:end)]);
 


return
