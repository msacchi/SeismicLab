function [r] = gauss_mixture(N,arg);
%GAUSS_MIXTURE: Compute a reflectivity using a Gaussian Mixture model.
%
% [r] = gauss_mixture(N,[sigma1,sigma2,p]);
%
% IN   N:      number of samples
%      arg=[sigma1,sigma2,p] where
%      sigma1: variance of distribution with probability p
%      sigma2: variance of distribution with probability 1-p
%      p:      mixing parameter (0,1)
%
% OUT  r:      reflectivity
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.
   

 sigma1 = arg(1);
 sigma2 = arg(2);
 p = arg(3);

 r1 = randn(N,1)*sqrt(sigma1);
 r2 = randn(N,1)*sqrt(sigma2);

 for k=1:N
   if rand(1,1) < p; r(k)=r1(k); else
                     r(k)=r2(k); end;
 end
