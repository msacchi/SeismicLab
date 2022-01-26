function [s] = stackgather(d,N);
%STACKGATHER: a program to stack (with normalization) one gather.
%             s(t) = [sum_x d(x,t) ]/N(t);
%
%  [s] = stackgather(d,N);
%
%  IN   d(nt,nh): data (gather)
%       N(nt,1):  normalization factor for each time 
%
%  OUT  s:        stack (a normalized spatial average of traces) 
%
%  If N is not provided, we divide by number of traces
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.




 [nt,nh] = size(d);

 M = nargin;

 if M<2; N = nh*ones(nt,1); end;

 s = sum(d,2);

 for k = 1:nt;

   if N(k,1)>0; s(k,1) = s(k,1)/N(k,1); else;
                s(k,1) = 0.; end;

 end;
