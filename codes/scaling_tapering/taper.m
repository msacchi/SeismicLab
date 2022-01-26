function [dout,win] = taper(d,perc_beg,perc_end);
%TAPER: Apply a triangular  taper to the beg/end of traces.
% 
%  IN   d(nt,nx):      data (columns are traces)
%       perc_beg:      percetage of data to be tapered at the beggining
%       perc_end:      "                                      end
%
%  OUT  dout(nt,nx):   data with beg/end tapered
%       win(nt):       the taper
%
%  Example
%
%    d = make_traces; wigb([d,taper(d,50,50)]);
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


 [nt,nx] = size(d);

 i1 = floor(nt*perc_beg/100)+1;
 i2 = floor(nt*perc_end/100)+1;
   
 win = [ [1:1:i1]/i1,ones(1,nt-i1-i2),[i2:-1:1]/i2];

 for k=1:nx
  dout(:,k) = d(:,k).*win';
 end;


