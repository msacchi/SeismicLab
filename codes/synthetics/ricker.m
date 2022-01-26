function [w,tw] = ricker(f,dt)
%RICKER: Ricker wavelet of central frequency f
%
%  [w,tw] = ricker(f,dt);
%
%  IN   f : central freq. in Hz (f <<1/(2dt) )
%       dt: sampling interval in sec  
%
%  OUT  w:  the Ricker wavelet
%       tw: axis
%
%  Example
%
%    [w,tw] = ricker(10,0.004);
%    plot(tw,w);
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


 nw=2.5/f/dt;
 nw=2*floor(nw/2)+1;
 nc=floor(nw/2);
 w = zeros(nw,1);

 k=[1:1:nw]';

 alpha = (nc-k+1).*f*dt*pi;
 beta=alpha.^2;
 w = (1.-beta.*2).*exp(-beta);

  if nargout>1;
    tw = -(nc+1-[1:1:nw])*dt;
  end
