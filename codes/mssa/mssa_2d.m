function [DATA_f,sing] = mssa_2d(DATA,dt,P,flow,fhigh,meth);
%MSSA_2D: 2D mssa filtering of a  2D gather  
%
%  [DATA_f] = mssa_2d(DATA,dt,P,flow,fhigh,meth);
%
%  IN   DATA:    The 2D array data, columns are traces DATA(t,x)
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
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.

        [nt,nx] = size(DATA);
        nf = 2*2^nextpow2(nt);

        DATA_FX_f = zeros(nf,nx);

 sing = [];
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

% Size of Hankel Matrix

        Lcol = floor(nx/2)+1;
        Lrow = nx-Lcol+1;

% Form level-1 block Hankel matrix 


        for j= ilow:ihigh

        M = zeros(Lrow,Lcol);

           for lc = 1:Lcol
             M(:,lc)  = DATA_FX_tmp(j,lc:lc+Lrow-1);
           end

% SVD deconcoposition with P largest singular values of Randomized SVD 

    if meth==1;     [U,S,V] = svds(M,10); Mout = U(:,1:P)*S(1:P,1:P)*V(:,1:P)'; end;
    if meth==2;                           Mout = rqrd(M,P); end;



% Sum along anti-diagonals to recover signal 


   Count = zeros(nx,1);
   tmp2 = zeros(nx,1);
   
     for ic = 1:Lcol;
      for ir = 1:Lrow;
       Count(ir+ic-1,1) = Count(ir+ic-1,1)+1;
       tmp2(ir+ic-1,1)  = tmp2(ir+ic-1,1) + Mout(ir,ic);
       end;
      end

       tmp2 = tmp2./Count;

     DATA_FX_f(j,:) = tmp2;

   end

% Honor symmetries

   for k=nf/2+2:nf
    DATA_FX_f(k,:) = conj(DATA_FX_f(nf-k+2,:));
   end

% Back to TX (the output)

   DATA_f = real(ifft(DATA_FX_f,[],1));
   DATA_f = DATA_f(1:nt,:);

return
