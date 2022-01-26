function [f,tf,x,tx,Med_Norm] = med(wbp,s,dt,Nf,mu,Updates,fun,dfun,arg);
%MED: Multichannel deconvolution via Minimum Entropy Deconvolution
%
% [f,tf,x,tx,Med_Norm] = med(wbp,s,dt,Nf,mu,Updates,fun,dfun,arg);
%
%  IN   wbp:      band-pass filter. For instance, [0.25,0.5,0.25]'
%       s:        the multichannel signal 
%       dt:       sampling interval in sec
%       Nf:       length of operator (length of the filter)
%       mu:       pre-whitening
%       Updates:  number of iterations
%       fun:      entropy function   (see Note below) 
%       dfun:     derivative of the entropy function  (see Note below) 
%
%  OUT  f:        filter 
%       tf:       time axis for filter 
%       x:        deconvolved signal
%       tx:       time axis for output deconvolved signal 
%       Med_Norm: norm 
%
%
%  Reference: Wiggins 1978, Minumum Entropy Deconvolution, Geoexploration, 16, 21-35. 
%
%             The Logartimic norm is in the following paper: 
%
%             Sacchi, Velis and Cominguez, 1995, Minimum entropy deconvolution with frequency‐domain constraints,
%             Geophysics, 59(6), 864-1017
%  Note: 
%
% for the logaritmic entropy norm use 
%
% [f,tf,x,tx,Med_Norm] = med(wbp,s,dt,Nf,mu,Updates,@funln,@dfunln,[]); 
%
% for varimax-style norm use: 
%
% [f,tf,x,tx,Med_Norm] = med(wbp,s,dt,Nf,mu,Updates,@funalpha,@dfunalpha,alpha); 
%
% if alpha =1 you have the Varimax norm and therefore, the original method proposed by Wiggins. 
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.

 [Nt,Nc] = size(s);

% Compute average Autocorrelation matrix

 Rs = zeros(Nf,Nf);

 for k=1:Nc;
  S = convmtx(s(:,k),Nf);
  Rs = Rs + S'*S;
 end

 Rs = Rs/Nc;

 mu = Rs(1,1)*mu/100;

% Initial MED Filter 

 f = zeros(Nf,1);
 f(floor(Nf/2),1) = 1;
 Nx = Nt+Nf-1; 

 for k=1:Nc
  x(:,k) = conv(s(:,k),f);
 end;

 Q = eye(Nf);

% Compute inverse and save it

 Matrix = Rs+mu*Q;
 Matrix = inv(Matrix);


 Med_Norm = [];

 for j = 1: Updates;
 
% Get RHS term

  vs = 0;
  g = zeros(Nf,1);
   for k=1:Nc;
    [b,v] = non_linear_output(x(:,k),fun,dfun,arg);
    S = convmtx(s(:,k),Nf);
    g = g + S'*b;
    vs = vs + v;
  end
  g = g/Nc;
  vs = vs/Nc;
  f = Matrix*g;
  Med_Norm(j) = vs;
 
% Update reflectivity

   x = conv2(s,f);

 end

% bandlimit output 
  x =  conv2(x,wbp,'same');
  x = x(1:Nt,:); 

 [x] = delay(s,x,Nf);
 tx = (0:1:Nt-1)*dt;
 tf = (0:1:Nf-1)*dt;

% Normalize ouput and filter 

 smax = max(s(:));
 xmax = max(x(:)); 
 c    = smax/xmax; 
 x    = c*x;
 f    = f/c; 

return
