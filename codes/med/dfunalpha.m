function g=dfunalpha(x,arg);
%dfunalpha: dericative of ALPHA norm (Wiggins' alpha=1) for MED deconvolution
%
alpha = arg;

eps= 1.e-10;
I = find(x>eps);

g = zeros(size(x));
g(I) = alpha.*x(I).^(alpha-1);
