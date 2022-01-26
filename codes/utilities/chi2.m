function Chi2 = chi2(dp,d,sigma);
%CHI2: Get Chi^2
%
%  [Chi2] = chi2(dp,d);
%
%  IN   dp:  estimated signal
%       d :  observed signal
%       sigma: noise level (std)  
%
%  OUT  Chi2:  Chi^2
%
%  Copyright (C) 2015, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


e = (dp-d).^2;

Chi2 = sum( e(:))/sigma^2; 

return;
