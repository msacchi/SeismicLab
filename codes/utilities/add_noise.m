function  [dout,sigma] = add_noise(din,SNR,L);
%ADD_NOISE: Add noise with a prescribed SNR.
%
% [dout] = add_noise(din,SNR,L);
% 
%  IN   din:     input data of any dimension
%       SNR:     desired SNR (SNR = Power of signal/Power of noise)
%       L:       lenght of Hamming window convolved 
%                in time to produce noise that low-pass. 
%                (L=1 is for white noise)
% 
%  OUT  dout:    data contaminated with noise 
%       sigma:   std of noise added to data 
%
%  Example:
%   
%      d = linear_events();
%      dn = add_noise(d,1,11); 
%      wigb([d,dn]); 
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


 N  = size(din);
 nt = N(1);
 X  = reshape(din,nt,prod(N(2:end)));
 h  = hamming(L);
 Noise = randn(size(X));
 Noise = conv2(Noise,h,'same'); 
 alpha = sqrt( sum(X(:).^2)/(SNR*sum(Noise(:).^2)) );
 noise_added = alpha*Noise;
 Y = X + noise_added;  
 dout = reshape(Y,[nt,N(2:end)]);

 sigma = std(noise_added(:)); 
 
return 
