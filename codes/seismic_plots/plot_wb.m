function  plot_wb(t,w,a);
%PLOT_WB: Horizontal plotting with bias.
%
%  IN    t: time axis
%        w: traces or wavelets in columns.
%        a: overlab in percentage
%
%
%  Example:
%            
%      [w1,tw] = rotated_wavelet(4./1000,1.,20.,0);
%      [w2,tw] = rotated_wavelet(4./1000,1.,20.,45);
%      [w3,tw] = rotated_wavelet(4./1000,1.,20.,90);
%      W =[w1,w2,w3];
%      subplot(221); plot_wb(tw,W,  0);
%      subplot(222); plot_wb(tw,W,50);
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


 w = w/max(max(abs(w)));
 [nt,nx] = size(w);
 x = 1:1:nx;
 plot(t,((a/100+1)*w)+ones(size(w))*diag(x),'b','LineWidth',1.);axis tight;
 return;
