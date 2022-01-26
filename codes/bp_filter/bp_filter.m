function  [o] =  bp_filter(d,dt,f1,f2,f3,f4);
%BP_FILTER: Apply a band-pass filter to a group of traces. The traces
%           can be given by a matrix or a  cube or a ND volume. Firt dimension is time. 
%
%  [o] = bp_filter(d,dt,f1,f2,f3,f4);
%
%  IN   d:    data (first dimention of volume is time)
%       dt:   sampling interval in sec
%       f1:   freq. in Hz
%       f2:   freq. in Hz
%       f3:   freq. in Hz
%       f4:   freq. in Hz
%
%   ^
%   |     ___________
%   |    /           \   Amplitude spectrum
%   |   /             \
%   |  /               \
%   |------------------------>
%      f1 f2        f3 f4
%
%  OUT  o:    output  (columns are traces)
%
%  Example: 
%
%    d=linear_events; 
%    dn = add_noise(d,0.5,3); 
%    dout = bp_filter(dn,0.004,1,3,30,40); 
%    wigb([d,dn,dout]);
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License. 


 Ndims = ndims(d);
 N     = size(d);
 nt    = N(1);
 x     = reshape(d,nt,prod(N(2:end)));
 
 k = nextpow2(nt);
 nf = 4*(2^k);

 i1 = floor(nf*f1*dt)+1;
 i2 = floor(nf*f2*dt)+1;
 i3 = floor(nf*f3*dt)+1;
 i4 = floor(nf*f4*dt)+1;

 up =  (1:1:(i2-i1))/(i2-i1);
 down = (i4-i3:-1:1)/(i4-i3);
 aux1 = [zeros(1,i1), up, ones(1,i3-i2), down, zeros(1,nf/2+1-i4) ];
 aux2 = fliplr(aux1(1,2:nf/2));

 c = 0; % zero phase (could apply rotations as well)
 F = ([aux1,aux2]');
 Phase = (pi/180.)*[0.,-c*ones(1,nf/2-1),0.,c*ones(1,nf/2-1)];
 Transfer = F.*exp(-i*Phase');


 X = fft(x,nf,1);

 for k = 1:prod(N(2:end))
  Y(:,k) = Transfer.*X(:,k);
 end

 o = ifft(Y,nf,1);

 o = real(o(1:nt,:));

 o = reshape(o, [nt, N(2:end)]);

 
return
