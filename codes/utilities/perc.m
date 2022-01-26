function x=perc(d,a);
%PERC: Get clip for imagesc
%
%  IN   d:       input data of any dimension
%       a:       clip (0,1)
% 
%  OUT  x:       the clip
%
%  Example:
%   
%      d = linear_events();
%      dn = add_noise(d,1,11); 
%      imagesc(dn,perc(dn,0.9)); 
%      colorbar;
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


x = a*[-1,1]*max(abs(d(:))); 

return


