function [e] = envelope(d);
%ENVELOPE: The envelope of a group of traces.
%
%  e = envelope(d);
%
%  IN   d(nt,np):  traces
%
%  OUT  e(nt,np):  envelopes
%
%
%  Example
%
%    [d,r,t] = make_traces; d=d(:,1); e = envelope(d); 
%    plot(t,d,t,e);
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.



[nt,np] = size(d);

  for k=1:np

  aux = d(:,k);
  u = hilbert(aux);
  e(:,k) = abs(u);

 end;



