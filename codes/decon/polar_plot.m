function polar_plot(z);
%POLAR_PLOT: Plot the roots of a wavelet in polar coordinates.
%
%  polar_plot(z)
%
%  IN   z: zeros of the wavelet, they can be computed
%          with zeros_wav.m
%
%  NOTE: some zeros migth end up outside the plot (check axis)
%
%  Example
%
%    wavelet = ricker(20,0.004); z = zeros_wav(wavelet); polar_plot(z);
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.

 
 x=cos(0:0.1:2*pi); y = sin(0:0.1:2*pi);
 
 II = find(abs(z)<1); z1=z(II);
 
 plot(x,y);hold on;plot(real(z),imag(z),'sk'); 
                  ;plot(real(z1),imag(z1),'*r'); 
axis equal;
 axis([-2,2,-2,2]);
return

