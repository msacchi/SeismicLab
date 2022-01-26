function [x] = laplace_mixture(N,arg);
%LAPLACE_MIXTURE: Series of random numbers computed from a Laplace mixture
%
% [x] = laplace_mixture(N,arg);
%
% IN   N:        Lenght of series 
%      arg:      arg = [lambda1,lambda2,p] where 
%      lambda1:  Laplace paramater for the first distribution
%      lambda2:  Laplace paramater for the second distribution
%      p:        Mixing parameter (0,1)
%
% OUT  x(N,1):   Series 
%
%
% Example from Walden and Hosken (Geophysical Prospecting, 1986):
%
% lambda1 = 0.007; lambda2 = 0.017; p = 0.24; 
% r = laplace_mixture(500, [lambda1, lambda2, p]);
% plot(r); title('Reflectivity');
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.





 lambda1 = arg(1);
 lambda2 = arg(2);
 p = arg(3);

 x1 = laplace(N,lambda1);
 x2 = laplace(N,lambda2);

 x = zeros(N,1);
for k=1:N
  if  p>rand(1,1); x(k,1) = x1(k);
         else;
                   x(k,1) = x2(k);
   end
end

%expected_var = (p)*2*lambda1^2 + (1.-p)*2*lambda2^2
%var(x) 

return


