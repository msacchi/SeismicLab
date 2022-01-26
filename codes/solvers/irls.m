function [x,J] = irls(x0,b,operator, Param, mu, max_iter_cgls, max_iter_irls, tol1,tol2)
%IRLS: Iterative Re-weighted Least-squares for sparse solutions. Find x that
%      minimize the l2-l1 cost function  J = ||A x- b||_2^2 + mu ||x||_1. 
%      The solution is estimated by solving a sequence of quadratic
%      problems via cglsw.m
%
%  [x,J] = irls(x0,b,operator, Param, mu, max_iter_cgls, max_iter_irls, tol1, tol2)
%
%  IN   x0:   Starting point  
%       b:    data 
%       A:    is given as an operator of the form 
%             out  = operator(in, Param, flag) 
%               flag =  1 ==> A
%               flag = -1 ==> A^T
%       mu:   trade-off parameter
%       max_iter_cgls:   maximum iterations for cgls
%       max_iter_irls:   maximum number of iterations
%       tol1:            tolerance for cgls (1.e-6)
%       tol2:            tolerance for IRLS outer loop (1.e-4)
%
%  OUT  x:    solution
%       J:    cost vesus iteration
%   
%  Notes:
%
%  tol1: Stoping criterion for cgls.m. It stops when the normalized 
%  l2 norm of the gradient of the quadratic cost function is 
%  less than tol1.
%
%  tol2: Outer loop stops  when the normalized l2  norm of the gradient  
%  of J is less thatn tol2.
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Modified by M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


 Wr = ones(size(b)); 
 Wx = ones(size(x0)); 
 u0 = zeros(size(x0)); 

 Diff = 99999.;

 for j = 1 : max_iter_irls

     u = cglsw(u0, b ,operator, Param, Wr, Wx, mu, max_iter_cgls, tol1, 0);
     x = Wx.*u;
     Wx = sqrt(abs(x) + 0.00001); 
     e = Wr.*(operator(x,Param,1)-b);
     
     J(j) = sum(abs(e(:)).^2)+mu*sum(abs(x(:)));

  if j>1; 
    Diff = abs(J(j) - J(j-1))/((J(j)+J(j-1))/2);
  end

  if Diff<tol2
   break
  end

 
 end

 fprintf(' IRLS ended after  %4.0f iterations of %4.0f  \n',  j,max_iter_irls)

return
  
