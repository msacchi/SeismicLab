clear

      N = [40,50] 
      nt = 340;     
      dx = [10.,10.]; 
      dt = 4./1000.;
      p  = [0.11, -0.4, 0.19,  0.3;
            0.2,  -0.1, 0.12,  0.3]
      t0 = [ 0.1, 0.4,  0.33, 0.5];
      A  = [-1.0, 1.0,  0.8, -0.9];
      f0 = 20.;
      opt = 'parabolic'

      [d] = data_cube(N,dt,f0,nt,dx,t0,p,A,opt);
      [d] = bp_filter(d,dt,2.,10.,20.,30);


      dn = add_noise(d,100,4); 
      P = 94;
      flow = 1
      fhigh = 31; 


     [nt,n1,n2] = size(d); 

     T = ones(size(dn));
     for i1=1:n1 
      for i2=1:n2;
       p=rand(1); if p<0.7; T(:,i1,i2) = 0; end
      end
     end
 

       dn = dn.*T; 
       [df] = mssa_3d(dn,dt,P,flow,fhigh,1);


% ---------------------

     [nt,n1,n2] = size(df); 

     Df = reshape(df,nt*n1,n2);
     Dn = reshape(dn,nt*n1,n2);

     [U,S,V] = svd(Df,'econ'); 
     P2 = 20;
     U = U(:,1:P2);

 % W (Df -  U A) 

     Dout = U*U'*Dn;

     d2 = reshape(Dout,nt,n1,n2);
     
    figure(1); clf;
    a1 = squeeze(dn(:,:,10)); 
    a2 = squeeze(df(:,:,10)); 
    a3 = squeeze(d2(:,:,10)); 

    wigb([a1,a2,a3]);
     


