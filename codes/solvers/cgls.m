function [x,J] = cgls(x0, b ,operator, Param, mu, max_iter, tol, prnt)
%CGLS: Solves for the minimum of J = || A x - b ||_2^2  + mu ||x||_2^2 
%      via the method of conjugate gradients for least-squares problems. The 
%      matrix A is given via an  operator  and apply on the flight by
%      user-defined function "operator" with parameters "Param".
%
%  x = cgls(x0, b ,operator, Param, mu, max_iter, tol, prnt)
%
%
%  IN   x0:       starting point  
%       A:        is given as an operator of the form 
%                 out = operator(in, Param, flag) 
%                  flag =  1 ==> A
%                  flag = -1 ==> A^T
%       b:        data 
%       mu:       trade-off parameter
%       max_iter: maximum number of iterations
%       tol:      stopping criterion (1.e-6)
%       prnt: 1:  print , 0: silent
%               
%  OUT  x:    solution
%       J:    cost versus iteration 
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Modified by M.D.Sacchi based on code by Saunders (https://web.stanford.edu/group/SOL/software/cgls/)
%
%  SeismicLab is licensed under the MIT License.

   x = x0;
   r = b - operator(x,Param,1);
   s = operator(r,Param,-1) - mu*x;
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

 
 while (k <max_iter) && (flag==0) ;
  
    q = operator(p,Param,1);
    delta = cgdot(q,q) + mu*cgdot(p,p);
    if delta == 0, delta = 1.e-10   ; end
    alpha = gamma/delta; 
    x = x + alpha*p;
    r = r - alpha*q;
    s = operator(r,Param,-1) - mu*x;

    gamma1  = cgdot(s,s);
     norms  = sqrt(gamma1);
       beta = gamma1/gamma;
      gamma = gamma1;
    
    p = s + beta*p;
    
       flag = (norms<=norms0 * tol);
       nres = norms / norms0;
          k = k+1;

      e = operator(x,Param,1)-b; 
      J(k) = sum( (abs(e(:))).^2  ) + mu*sum( (abs(x(:))).^2 );

 if prnt, fprintf(form, k, norms, nres); end
 
 end; 

 
 % Diagnostics
 
 if  k == max_iter;   flag = 2; end;
 if prnt, 
     
    
            fprintf( ' ============================================== \n');
if flag==1; fprintf( ' ====== CGLS converged before max_iter ======== \n'); end
if flag==2; fprintf( ' ====== CGLS reached max_iter ================= \n'); end
            fprintf( ' ============================================== \n \n \n');
 end
 
 fprintf(' CGLS ended after  %4.0f iterations of %4.0f  \n',  k,max_iter)
 return;
  
