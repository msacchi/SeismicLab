function [x] = cglsw(x0, b ,operator, Param, Wr, Wx, mu, max_iter, tol, prnt)
%cgls:  Solves for the minimum of  J = ||Wr ( A.Wx.x - b) ||_2^2  + mu ||x||_2^2 via  
%       the method of conjugate gradients for least-squares problems. The
%       matrix A is given in operator form. Both A and A^T are  applied on the
%       the flight via a user-defined function "operator" with parameters "Param".
%
%  x = cglsw(x0, b ,operator, Param, Wr, Wx,  mu, max_iter, tol, prnt)
%
%
%  IN   x0:         Starting point
%       A:          is given as an operator of the form
%       Wr:         element-to-element weights for residuals 
%       Wx:         element-to-element weights for model parameters
%
%                   out=operator(in, Param, flag)
%                   flag =  1 ==> A
%                   flag = -1 ==> A^T
%
%       b:          data
%       mu:         trade-off parameter
%       max_iter:   maximum number of iterations
%       tol:        stopping criterion (1.e-6)
%       prnt:       prnt=1 print diagnostics
%                   prnt=0 silent
%               
%  OUT  x:          solution
%
%  Note: This program is needed by irls.m
%    
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Modified by M.D.Sacchi based on the code by Saunders (https://web.stanford.edu/group/SOL/software/cgls/)
%
%  SeismicLab is licensed under the MIT License.

 b = Wr.*b;

   x = x0;
   r   = b - Wr.*operator(Wx.*x,Param,1);
   s   = Wx.*operator(Wr.*r,Param,-1) - mu*x;
   p = s;

   gamma  = cgdot(s,s);
   norms0 = sqrt(gamma);           % norm of the gradient is used to stop 
   k      = 0;
   flag   = 0;

if prnt 
    
    fprintf( ' ============================================== \n');
    fprintf( ' ================= CGLS ======================= \n');
    fprintf( ' ============================================== \n');
    
    head = '     k           |grad|       |grad|/|grad_0|        '; 
    form = '   %3.0f       %12.5g           %8.3g              \n';
    disp('  ');   disp(head);
    fprintf(form, k, norms0,1)

end

 
 while (k <max_iter) && (flag==0);
  
    q = Wr.*operator(Wx.*p,Param,1);
    delta = cgdot(q,q) + mu*cgdot(p,p);
    if delta == 0, delta = 1.e-10; end
    alpha = gamma/delta; 
    x = x + alpha*p;
    r = r - alpha*q;
    s = Wx.*operator(Wr.*r,Param,-1) - mu*x;

    gamma1 = cgdot(s,s);
    norms = sqrt(gamma1);
    beta = gamma1/gamma;
    gamma = gamma1;
    
    p = s + beta*p;
    
    flag = (norms<=norms0 * tol);
    nres = norms / norms0;

    k = k+1;

 if prnt, fprintf(form, k, norms, nres); end
 
 end; 

 
 % Diagnostics
 
 if  k == max_iter; flag = 2; end

if prnt,
     

               fprintf( ' ============================================== \n');
   if flag==1; fprintf( ' ====== CGLS converged before max_iter ======== \n'); end
   if flag==2; fprintf( ' ====== CGLS reached max_iter ================= \n'); end
               fprintf( ' ============================================== \n \n \n');

end

return
  
