function g=funln(x,arg);
%FUNLN: function for logarithmic MED deconvolution
p = arg;
g = zeros(length(x),1);
I = find(abs(x)>1.e-10);
g(I) = log(x(I));
