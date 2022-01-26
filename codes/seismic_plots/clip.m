function  dc  = clip(d,cmin,cmax);
%CLIP: A program to clip seismic traces
%
%  dc  = clip(d, cmin, cmax);
%
%  IN   d:          data to be clipped, one trace, many traces in a matrix, cube so one
%                   first dimension is time
%       cmin:       lower clip in % (90% means clip at the 90% negative amplitude)
%       cmax:       lower clip in % (90% means clip at the 90% positive amplitude)
%
%  OUT  Dc:         data after being clipped
%
%  Example: 
%
%  d = sin(2*pi*.02*[1:1:200]); plot(clip(d,90,90));
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


if nargin==3;
dmax = (cmax/100)*max(d(:));
dmin = (cmin/100)*min(d(:));
I = find(d>dmax); d(I)=dmax;
J = find(d<dmin); d(J)=dmin;
end
if nargin==2;
dmax = (cmin/100)*max(d(:));
dmin = -dmax;
I = find(d>dmax); d(I)=dmax;
J = find(d<dmin); d(J)=dmin;
end

dc = d;
  return;
