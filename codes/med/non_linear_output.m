function [b,h] = non_linear_output(y,fun,dfun,arg);
%NON_LINEAR_OUTPUT:evaluation of the non-linear desired output
%                  for MED. 

  N = length(y);
aux = y.^2;
 c1 = sum(aux)/N;
  q = aux/c1;

  F =  feval(fun,q,arg);
 dF = feval(dfun,q,arg);
  G = F+q.*dF;
 c2 = sum(G.*q)/N;

  b = G.*y/c2;
 c3 = N*feval(fun,N,arg);
  h = sum(q.*F)/c3;;
