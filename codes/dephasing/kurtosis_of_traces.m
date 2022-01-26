function [K] = kurtosis_of_traces(x)
%KURTOSIS: Kurtosis of a one or more time series.
%
%  [K] = kurtosis_of_traces(x);
%
%  IN   x:      data (vector or matrix)
%
%  OUT  K:      Kurtosis
%
%
%  Kurtosis is defined as K = E(x^4)/ (E(x^2))^2;
%
%    K = 3 for a Gaussian series 
%
%  There is a different definiton k' = K-3
%
%    k' = 0 for a Gaussian series
%
%  k' is called the Kurtosis Excess. This matlab 
%  function computes K not k'.
%
%  Reference: Longbottom, J., Walden, A.T. and White, R.E. (1988) Principles and 
%             application of maximum kurtosis phase estimation. Geophysical 
%             Prospecting, 36, 115-138.
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


 N = ndims(x);

[n1,n2]=size(x);
 Ns = n1*n2;

 if N==1; 
  sum1 = sum(x.^4)/Ns;
  sum2 = (sum(x.^2)/Ns)^2;
 end;
 
 if N==2; 
  sum1 = sum(sum(x.^4)/Ns);
  sum2 = (sum(sum(x.^2))/Ns)^2;
 end;

 K = sum1/sum2;

 return;
