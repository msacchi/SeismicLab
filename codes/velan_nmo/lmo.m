function [dout] = lmo(d,dt,h,v);
%LMO: A program for linear moveout
%
% [dout,M] = nmo(d,dt,h,v);
%
%  IN   d(nt,nh):  data (gather)
%       dt:        sampling interval in secs
%       h(nh):     vector of offsets in meters
%       v:         linear moveout velocity   
%
%  OUT  dout:      data after LMO correction
%
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.

 

% Intepolate t0,v pairs 

  [nt,nh] = size(d);

  dout = zeros(size(d));

  for ih = 1:nh;

  arg =   h(ih)/v;

  time = arg;

  its = time/dt+1;
  it1 = floor(time/dt+1);
  it2 = it1+1;
  a = its-it1;

   if it2 <= nt; dout(it,ih) = (1-a)*d(it1,ih)+a*d(it2,ih); end;

  end

 return;
    
