function sgray(alpha)
%SGRAY: Non-linear transformation of a gray colormap.
%       Similar to clipping.
%
%  sgray(alpha)
%
%  IN  alpha: degree of BW color scale clustering (try 0.5)
%
%  Example: 
%
%  d = linear_events; imagesc(d); sgray(0.5);
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.




i0=32;
i = 1:64;
t = (atan((i-i0)/alpha))';
s = t(64);
t = (t - min(t))*1./(max(t) -min(t));


m(1:64,2) = 1-t;
m(:,1) = 1-t;
m(:,3) =1-t;
colormap(m);
