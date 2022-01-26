function [x] = laplace(N,lambda);
%LAPLACE: Random series with Laplace distribution.
%         Distribution  f(x) = 1/2/lambda * e{-|x|/lambda}}
%
% [x] = laplace(N,lambda);
%
% IN   N:      Lenght of series 
%      lambda: Laplace parameter 
%
% OUT  x(N,1): Series 
%
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.
  



 y = rand(N,1);

 for k=1:N;
   if y(k)<=0.5; x(k)= lambda*log(2*y(k)     );    else
                 x(k)=-lambda*log(2*(1-y(k)));
  end;
 end;

 return;
