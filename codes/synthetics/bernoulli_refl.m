function r = bernoulli_refl(N,arg);
%BERNOULLI: Random numbers with Bernoulli distribution
%
% [r] = bernoulli(N,lambda,sigma);
%
% IN   N:    length of series 
%      arg = [lambda,sigma]); 
%      lambda: Occurrence of a non-zero sample (0,1) 
%      sigma: standard error for non-zero samples 
%
% OUT  r(N,1): random numbers (series) 
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.

 lambda = arg(1); 
  sigma = arg(2); 

 r = zeros(N,1);

 for k=1:N
  if rand>lambda; r(k,1)= sigma*randn; end
 end;




