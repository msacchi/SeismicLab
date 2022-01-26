function pimage(x,z,d);
%PIMAGE: High quality image for ppt presentations etc etc. 
%        Go to "file" "export setup" to export with proper labels, 
%        line width, size and other attributes.
%
%  IN  x: x-axis, x(nx)
%      z: z-axis, z(nz)
%      d: data, d(nz,nx)
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.



 pcolor(x,z,d);
 shading interp;
 axis ij; 


 set(gca,'ydir','reverse','xaxislocation','top','yaxislocation','left','layer','top','linewidth',2,'tickDir','out','box','on')

 % Move the plot down to make space at the top

  pos=get(gca,'position');
  set(gca,'position',[pos(1) 0.05 pos(3) pos(4)])

 colormap(seismic(1));

 return;


