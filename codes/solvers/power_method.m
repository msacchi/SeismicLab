function value = power_method(x0,Hop,PARAM);
%POWER_METHOD: Power iteration method to computes max eigenvalue of H'H 
%              where H and H' are given by Hop. This function is 
%              needed to evaluate the step parameter of FISTA.
%
%  IN   x0:    initial seed with dimensions such that H'*H*x0 does not abort 
%       Hop:   liner operator that encapsulates H and H' such that 
%              Hop(x,PARAM, 1) = H x
%              Hop(y,PARAM,-1) = H'y 
%       PARAM: set of parameters needed by Hop, PARAM is a structure 
%  
%  OUT  value: maximum eigenvalue of H'H
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.%

 x = x0;

 for k = 1:10;
    
    aux = Hop(x,PARAM,1);                  % aux = L x
    y = Hop(aux,PARAM,-1);                 % y = L'aux = L'L x 
    n = norm(y(:));                         
    x = y/n;                               % x = y/|x|
    value = n;

     fprintf('%6.0f %10.4f\n',k-1,value);

 end;

value = n;
