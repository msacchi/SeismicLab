function [dout] = fxy_eigen_images(d,dt,P,flow,fhigh,meth);
% EIGEN_IMAGED_3D: Denoising in the FXY using Trickett's eigenimage filtering
%
%  [dout] = fxy_eigen_images(d,dt,P,flow,fhigh,meth);
%
%  IN   d:       The 3D array data, columns are traces d(t,x,y)
%       dt:      sampling interval in sec
%       P:       size of subspace for noise attenuation (desired rank approx number of dips)
%       flow:    min  freq. in the data in Hz
%       fhigh:   max  freq. in the data in Hz
%       meth:    1 = standart SVD, 2 = Randomized SVD, 3= CUR
%
%  OUT  dout:    filtered data
%
%  Example:
%     
%        d = data_cube();
%        dn = add_noise(d,1.,3);
%        df = fxy_eigen_images(dn,0.004,4,1.0,120.0,1);  
%        wigb([d(:,:,3), dn(:,:,3), df(:,:,3)],2); 
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


 [nt,nx,ny] = size(d);
 nf = 2*2^nextpow2(nt);

 dout = zeros(nt,nx,ny);
 Dout = zeros(nf,nx,ny);

% First and last sample of the DFT.

 ilow  = max(floor(flow*dt*nf)+1,1);
 ihigh = min(floor(nf/2)+1,nf/2+1);

% Transform to FX

D = fft(d,nf,1);

for j= ilow:ihigh

 M = squeeze(D(j,:,:));

 if meth==1;     [U,S,V] = svds(M,P); Mout = U*S*V'; end;
 if meth==2;                          Mout = rand_svd(M,P,2*P,2); end; 
 if meth==3;                          Mout = cur(M,P); end; 

 Dout(j,:,:) =  Mout(:,:); 

end;

% Honor symmetries

 for k=nf/2+2:nf
  Dout(k,:,:) = conj(Dout(nf-k+2,:,:));
 end

% Back to TX to obtain the output

 dout = (ifft(Dout,[],1));
 dout = dout(1:nt,:,:);

return
