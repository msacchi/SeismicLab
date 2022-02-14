function [DATA_f] = mssa_3d(DATA,dt,P,flow,fhigh,meth);
%MSSA_3D: applies mssa filtering to a 3D cube. 
%
%  [DATA_f] = mssa_3d(DATA,dt,P,flow,fhigh,meth);
%
%  IN   DATA:    The 3D array data, columns are traces DATA(t,x,y)
%       dt:      sampling interval in sec
%       P:       size of subspace for noise attenuation (desired rank is the number of dips)
%       flow:    min  freq. in the data in Hz
%       fhigh:   max  freq. in the data in Hz
%       meth:    1 = standart SVD, 2 = Randomized SVD
%
%  OUT  DATA_f: filtered data
%
%  Reference: Oropeza and Sacchi, 2011, Simultaneous seismic de-noising and reconstruction via Multichannel
%             Singular Spectrum Analysis (MSSA), Geophysics 76(3)
% 
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.

        [nt,nx,ny] = size(DATA);
        nf = 2*2^nextpow2(nt);

        DATA_FX_f = zeros(nf,nx,ny);

% First and last samples of the DFT.

        ilow  = floor(flow*dt*nf)+1;

        if ilow<1;
            ilow=1;
        end;

        ihigh = floor(fhigh*dt*nf)+1;

        if ihigh > floor(nf/2)+1;
            ihigh=floor(nf/2)+1;
        end

% Transform to FX

        DATA_FX_tmp = fft(DATA,nf,1);
        samp = 1:1:nf;
        f = (samp-1)/(nf*dt);
        f = f(1:ihigh);

% Size of the Hankel Matrix for y.

        Ncol = floor(ny/2)+1;
        Nrow = ny-Ncol+1;

% Size of Hankel Matrix of Hankel Matrixes in x.

        Lcol = floor(nx/2)+1;
        Lrow = nx-Lcol+1;


% Form level-2 block Hankel matrix 

        for j= ilow:ihigh

          M = 0;
           for lc = 1:Lcol
            for lr = 1:Lrow

             tmp_fx(lr+lc-1,:)  = squeeze(DATA_FX_tmp(j,lr+lc-1,:)).';

              for ic = 1:Ncol;
                for ir = 1:Nrow;

                  M(((lr*Nrow)-Nrow+ir),((lc*Ncol)-Ncol+ic)) = tmp_fx(lr+lc-1,ir+ic-1);

                end
              end
            end
           end

   if j==ilow;
     fprintf(' ---- Size of Block Hankel Matrix:   %d x  %d \n',size(M));
   end

% SVD deconcoposition with P largest singular values of Randomized SVD 

    if meth==1;     [U,S,V] = svds(M,P); Mout = U*S*V'; 
    if  j==10;  s10 = diag(S); save('sing_10','s10'); end
    if  j==20;  s20 = diag(S); save('sing_20','s20'); end
    if  j==30;  s30 = diag(S); save('sing_30','s30'); end
    if  j==40;  s40 = diag(S); save('sing_40','s40'); end
    if  j==50;  s50 = diag(S); save('sing_50','s50'); end
    if  j==60;  s60 = diag(S); save('sing_60','s60'); end
    end;
    if meth==2;                          Mout = rqrd(M,P); end;

% Sum along anti-diagonals to recover signal 

   Count = zeros(ny,nx);
     tmp = zeros(ny,nx);
    tmp2 = zeros(ny,nx);

    for lc = 1:Lcol
     for lr = 1:Lrow
      for ic = 1:Ncol;
       for ir = 1:Nrow;

            Count(ir+ic-1,lr+lc-1) = Count(ir+ic-1,lr+lc-1)+1;
            tmp(ir+ic-1,lr+lc-1)  = tmp(ir+ic-1,lr+lc-1) + Mout(((lr*Nrow)-Nrow+ir),((lc*Ncol)-Ncol+ic));

       end;
      end

            tmp2(:,lr+lc-1) = tmp(:,lr+lc-1)./Count(:,lr+lc-1);
            DATA_FX_f(j,lr+lc-1,:) = tmp2(:,lr+lc-1).';

     end
    end
   end

% Honor symmetries

   for k=nf/2+2:nf
            DATA_FX_f(k,:,:) = conj(DATA_FX_f(nf-k+2,:,:));
   end

% Back to TX (the output)

   DATA_f = real(ifft(DATA_FX_f,[],1));
   DATA_f = DATA_f(1:nt,:,:);

return
