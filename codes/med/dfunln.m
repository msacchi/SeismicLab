function g=dfunln(x,arg);
%DFUNLN: derivative of the log norm for MED deconvolution
%
%
p = arg;
g = zeros(length(x),1);
I = find(abs(x)>1.e-10);
g(I) = 1./x(I);
