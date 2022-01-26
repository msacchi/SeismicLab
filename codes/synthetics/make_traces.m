function [s,r,t] = make_traces(Nr,Ntraces,dt,w,r_fun,arg,rho);
%MAKE_TRACES: Makes an emsemble of traces.
%
%  [s,r,t] = make_traces(Nr,Ntraces,dt,w,r_fun,arg,rho);
%
%  IN   Nr:      mumber of samples of the reflectivity
%       Ntraces: number of traces
%       dt:      sampling interval in secs
%       w:       wavelet  (a column)
%       arg:     argument vector of r_fun 
%       rho:     degree of similarity from trace to trace 
%                  (rho=0 - no correlation between traces)
%                  (rho=1 first trace is reapeted Ntraces times)
%
%       r_fun options are 'bernoulli', 'laplace_mixture', 'gauss_mixture', 
%
%  OUT  s:       seismic traces (lenght = Nr+Length of wavelet-1)
%       r:       reflectivity sequences (lenght = Nr)
%       t:       time axis of r and s
%
%  Example
%
%    w = ricker(40,0.002);
%    [s,r,t] = make_traces(100,20,0.002,w,'gauss_mixture',[1.,0.01,0.1],0.1);
%    wigb(s);
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


if nargin==0;

 f0 = 25;
 dt = 2./1000;
 w = ricker(f0,dt);
 Nr = 200;
 Ntraces = 20;
 r_fun = 'gauss_mixture';
 arg = [0.0001,0.1,0.9];
 rho = 0.82;

end;
 
 nw = length(w);

  r = zeros(Nr,Ntraces);
 rr = feval(r_fun,Nr,arg);

 for k = 1:Ntraces
  r(:,k) = rho*rr + (1-rho)*feval(r_fun,Nr,arg);
 end;
 
 ss = conv2(r,w);
 n2 = floor(nw/2)+1;
 s = ss(n2:n2+Nr-1,:);

 Nt = length(s);
 if nargout>1;
  t = (0:1:Nt-1)*dt;
 end

s = taper(s,10,10);

