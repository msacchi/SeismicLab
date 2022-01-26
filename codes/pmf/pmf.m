function [D] = pmf(T,Dobs,a,niter,p,meth,ty);
% pmf: Paralell matrix factorization for tensor completion 
%      modified by MDS. Notice that  I used the random QR decomposition
%      to perform rank reduction in each mode rather than the matrix factorization
%      descrived in 
%
%      J Gao, A Stanton, MD Sacchi, 2015,
%      Parallel matrix factorization algorithm and its application to 
%      5D seismic reconstruction and denoising:
%      Geophysics, 80 (6)
%  
%  D = pmf(T,Dobs,a,niter,p);
%
%  IN     T:       sampling operator
%         Dobs:    data with missing entries (3rd or 4th order tensor)
%         a:       tradeoff parateter (0.8-0.95)
%         niter:   max iterations
%         p:       multi-rank that is applied to each  unfolding
%                  p(3) or p(4)  depending on order of tensor
%         method:  1 = Randomized QR decomposition
%         method:  2 = Alternating Least-squares factorization
%         ty:      1 = Gaussian noise
%         ty:      2 = Non-Gaussian (Erratic noise) case
%
%  OUT    D:  tensor after completion
%
%  Copyright (C) 2016, Signal Analysis and Imaging Group
%  Author: Mauricio D Sacchi
%


 D = Dobs;  n = size(D);
 N_modes = ndims(Dobs);
 Ones = ones(size(Dobs));
 A = ones(size(Dobs));
 epsi = 0.001;

 for k = 1:niter;
    
    C = zeros(size(Dobs));
    A = zeros( size(C)); 
    
 if meth==1; 
     for j = 1:N_modes;
        Y = doit_1(D,p,n,j);     % ramdomized-QR  
        C  = C + Y;
    end;
 else

    for j = 1:N_modes;
        Y = doit_2(D,p,n,j);    % Alternating LS
        C  = C + Y;
    end;
 end
if ty==1;      D = (Ones-a.*T).*C/4 + a.*Dobs; end
if ty==2; 
    D = (Ones-A.*T).*C/4 + A.*Dobs;
    E = T.*(Dobs-D);
    A = 1./(1.+4*a*sqrt(abs(E).^2+epsi^2));
end;
end;

return

  
function  C = doit_1(D,p,n,j); 
%
% Performs rank reduction of unfolded tensor in mode j via 
% randomized qr decomposition 

   A = unfold(D,j);
   X = rqrd(A.',p(j));      % this part is different to Gao's paper
   C = fold(X.',j,n);       % ask me about .' (it is make rqrd go fast !!)

return


function  C = doit_2(D,p,n,j); 
%
% Performs rank reduction of unfolded tensor in mode j via
% alternating least-sqaures factorization 
% 
% Only 3 iterations of ALS are run 

   M = unfold(D,j);
   n2 = size(M,2);
   B = randn(p(j),n2);

  for k = 1:3;
   A = M*B';
   B = (A'*A+0.00001*eye(p(j)))\(A'*M); 
  end

 C = fold(A*B, j,n);

return
