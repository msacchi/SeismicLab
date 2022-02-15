function [dout] = completion(T,d,a,f1,f2,dt,niter,p,meth,ty);
% completion: Engine for tensor completion in f-x domain.
%             Works for 3rd and 4th order tensors.
%
%
%  IN     T:       sampling operator in space T(:,:,:) or T(:,:,:,:)
%         d:       data can be d(:,:,:,:)  first dimension is time
%                  or d(:,:,:,:,:) first dimension is time. In other
%                  words, the input is a 5D t-x or a 4D t-x tensor
%         a:       Trade-off paramter (a=0.8-0.95)
%         dt:      sampling interval in seconds
%         f1:      first freq. to process in Hz
%         f2:      last freq. to process in Hz
%         niter:   number of iterations
%         p:       multi-rank p(3) or p(4)
%         meth:    1: rQR decompositon
%                  2: iterative alternative LS solution (original paper)
%         ty:      1: Gaussian noise
%         ty:      2: Non-Gaussian (Erratic) noise 
%
%  OUT    X:       Reconstructed tensor  
%
%  Note:  d(t,n1,n2,n3,n4) is a 5D seismic data tensor but
%         in reality we will do 4D completion because the tensor
%         is first map to f-x and one frequency at the time
%         is used in the 4D completion process.
%         * if d(t,n1,n2,n3) then we have 3rd order tensor 
%         completion in space
%
%  Copyright (C) 2013, Signal Analysis and Imaging Group
%  Author: N Kreimer & Mauricio D Sacchi
% 

 M = ndims(d);    % Number of dimensions
 N = M-1;         % Number of spatial dimensions 

 vect = size(d);

 nt = vect(1);
 nf = 2*2^nextpow2(nt);


% Define min-max frequency indeces

 if1 = floor(f1*nf*dt)+1;
 if2 = floor(f2*nf*dt)+1;

% Transform from t-x to f-x


 D = fft(d,nf,1);
 Dout = zeros(size(D));

    fprintf('============================================ \n');
    fprintf('Start N-D reconstruction via PMF  \n')
    fprintf('Total number of grid points:        %3i \n', prod(size(T)));
    fprintf('Total number of empty grid points:  %3i \n', sum(T(:)));
for ifreq = if1:if2;

if mod(ifreq,10)==0;    fprintf('Frequency number %3i of %3i \n',ifreq,if2-if1+1);end
    
    if (N==4) aux = squeeze(D(ifreq,:,:,:,:)); end;
    if (N==3) aux = squeeze(D(ifreq,:,:,:));   end;
    
    [aux] = pmf(T,aux,a,niter,p,meth,ty);
    
    if (N==4); Dout(ifreq,:,:,:,:)      = aux;
        Dout(nf-ifreq+2,:,:,:,:) = conj(aux); 
    end;
    
    if (N==3); Dout(ifreq,:,:,:)        = aux;
        Dout(nf-ifreq+2,:,:,:) = conj(aux); 
    end;
    
end;

 dout = ifft(Dout,nf,1);

 if N==4; dout = real(dout(1:nt,:,:,:,:)); end
 if N==3; dout = real(dout(1:nt,:,:,:));   end

return;
