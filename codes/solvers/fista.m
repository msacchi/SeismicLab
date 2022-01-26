function [x,J] = fista(x0,y,Hop,PARAM,mu,Nit,tol)
%FISTA: Solves the l2-l1 problem via Fast Iterative Shrinkage-Thresholdng Algorithm 
%       Given a linear operator H and it's adjoint H', the algorithm minimizes
%       J = ||H x - y||_2^2 + mu ||x||_1, where H is the  linear operator
%       encapsulated in Hop
%
%  [x, J] = fista(x0, y, Hop, PARAM, mu,  Nit, tol)
%
%  IN.    x0:     initial sol (zeros) to get size of x0
%         y:      data
%         Hop:    linear operator
%         PARAM:  parameters to run Hop and it's adjoint
%         mu:     trade-off parameter
%         x0:     initial sol just to get size of x
% 
%  Out.   x:      sparse solution
%         J:      cost functon vs iteration
%
%
%  Reference:  Beck and Teboulle,  2009, A Fast Iterative Shrinkage-Thresholding Algorithm
%              for Linear Inverse Problems∗ SIAM J. Imaging Science,  Vol 2 (1), 183-202
% 
%  Note: The function is called as follows 
%
%         [x,J] = fista(x0, y, @foo, PARAM, mu, Nit, tol)
%
%        where foo.m is a function 
%
%         out = foo(in, PARAM, 1)  if forward operator and equivalent to out = H* in
%         out = foo(in, PARAM,-1)  if adjoint operator and equivalent to out = H'*in 
%
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Modified by M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


% Use power method to get step lenght alpha 

 x0 = randn(size(x0)); 
 alpha  = 1.05*power_method(x0,Hop, PARAM);

 
 J = zeros(1, Nit);                  % Objective function
 x = zeros(size(x0));
 T = mu/(2*alpha);
 
 t = 1;                           
 yk = x;

 Diff = 10000;
 for k = 1:Nit

    tmpx = x;                       
    Hx = Hop(yk,PARAM,1);
    x = thresholding(yk + (Hop( y-Hx,PARAM,-1))/alpha,'s', T);
    J(k) = sum(abs(Hx(:)-y(:)).^2) + mu*sum(abs(x(:)));


    if k>1; Diff = abs(J(k) - J(k-1))/((J(k)+J(k-1))/2);
 end

    if Diff<tol
             break
    end
  
   % FISTA acceleration part

    tmpt = t;
    t = (1+sqrt(1+4*t^2))/2;
    yk = x + (tmpt-1)/t*(x-tmpx);   
    
 end

 fprintf('FISTA ended after  %4.0f iterations of %4.0f  \n',  k,Nit) 
 
return 
