function [f,o] = predictive(w,NF,L,mu);  
%PREDICTIVE: Predictive deconvolution filter.
%
%  [f,o] = predictive(w,NF,L,mu)
%
%  IN   w:  the wavelet or input trace 
%       NF: lenght of the inverse filter
%       L:  Prediction distance
%       mu: Percentage of prewhitening    
%
%  OUT  f:  the filter
%       o:  the ouput or convolution of the filter with 
%           the wavelet or trace 
%
%  Example:
%
%     w = ricker(30,0.004);
%     p = zeros(100,1); p(20,1)=4; p(40,1)=-3.; p(60,1)=2.; p(80,1)=-1.;
%     trace = conv(p,w);
%     [f,o] = predictive(trace,35,20,0.1);
%     subplot(211); plot(trace); title('Trace')
%     subplot(212); plot(o); title('Trace after predictive decon')
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


NW = max(size(w));           % Lenght of the wavelet

NO = NW+NF-1;                % Leght of the output      

[mc,mr]=size(w);
if mc <= mr; w = w'; end;

b = zeros(NO,1);
b = [w(L+1:NW)',zeros(1,NO-NW+L)]';

C = convmtx(w,NF);           % Convolution matrix 
 
R = C'*C;                    %  Toeplitz Matrix
r0 = R(1,1);
R = R + r0*mu/100;
rhs = C'*b;                  %  Right hand side vector 
f = R\rhs;                   %  Filter 

if nargout == 2
if L==1; f = [1,-f']'; else
f = [1,zeros(1,L-1),-f']';
end;

o = conv(f,w);               %  Actual output
o = o(1:length(w));
end


return

