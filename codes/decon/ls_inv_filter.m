function [f,o] = ls_inv_filter(w,NF,Lag,mu);  
%LS_INV_FILTER: Least-squares inverse filter of a wavelet.
%
%  [f,o] = ls_inv_filter(w,NF,Lag,mu)
%
%  IN   w:   the wavelet
%       NF:  lenght of the inverse filter
%       Lag: the position of the spike in the desired output     
%            Lag=1 for minimum phase wavelets 
%       mu:  prewhitening as s percentage of the zero lag autocorrelation
%            coefficient of the wavelet
%
%  OUT  f:   the filter
%       o:   the ouput or convolution of the filter with 
%            the wavelet 
%
%  Example:
%  
%    w = [4;-3;1;0;0;0];
%    [f,o] = ls_inv_filter(w,10,1,0.1);
%    subplot(331); stem(w); title('wavelet');
%    subplot(332); stem(f); title('filter');
%    subplot(333); stem(o); title('conv of filter with wavelet')
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.

NW = length(w);              % lenght of the wavelet

NO = NW+NF-1;                % Leght of the output      

[mc,mr]=size(w);
if mc <= mr; w = w'; end;

b = [zeros(1,NO)]';          % Desire output 
b(Lag,1) = 1.;               % Position of the spike

C = convmtx(w,NF);           % Convolution matrix 
 
R = C'*C;                    % Toeplitz Matrix
r0 = R(1,1);
R = R + r0*mu/100;

rhs = C'*b;                  % Right hand side vector 
f = R\rhs;                   % Filter 
                             % Inversion should have been done with 
                             % a fast solver (Levinson recursion) 

if nargout == 2
o = conv(f,w);               %  Actual output
end

return

