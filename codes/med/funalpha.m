function g=funalpha(x,arg);
%funalpha: ALPHA norm (Wiggins' alpha=1) for MED deconvolution
alpha = arg;
g = x.^alpha;
