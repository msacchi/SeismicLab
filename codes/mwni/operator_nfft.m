function [output] = operator_nfft(input,T,NDim,N,K,option);
% OPERATOR_NFFT: a wrapping for the Forward and Adjoint ND Fourier operator implemented via FFTs
% 
% References: 
%
%      Minimum weighted norm interpolation of seismic records
%      Bin Liu, Mauricio D. Sacchi
%      GEOPHYSICS Nov 2004, Vol. 69, No. 6, pp. 1560-1568
%
% IMPORTANT: The flag option dictates what the operator will do.
%
% option = 1   --> This is the forward operator also called synthesis operator. This operator synthesis
%                  spatial data from the wavenumber domain. In other words, it is an inverse FFT
%                  Therefore input = data in wavenumber
%                            ouput = data in space  
% option =-1   --> This is the adjoint operator also called the analysis operator. This operator
%                  takes data in space and map them to wavenumber  
%                  Therefore input = data in space
%                            ouput = data in wavenumber
%
%      T is the sampling operator with 1 when the trace is alive and 0 when the trace is dead 
%      T(n1)                                if NDim = 2
%      T(n1,n2)                             if NDim = 3
%      T(n1,n2,n3)                          if NDim = 4
%      T(n1,n2,n3,n4)                       if NDim = 5
%      
%      NDim: number of dimensions of data_binned including time (2,3,4,5)
%
%      N:    number of samples in each spatial dimension
%      N = n1              if NDim = 2
%      N = [n1,n2]         if NDim = 3
%      N = [n1,n2,n3]      if NDim = 4
%      N = [n1,n2,n3,n4]   if NDim = 5
%
%      K:    number of wavenumbers in each direction 
%      K = k1              if NDim = 2
%      K = [k1,k2]         if NDim = 3
%      K = [k1,k2,k3]      if NDim = 4
%      K = [k1,k2,k3,k4]   if NDim = 5
%
% IMPORTANT: The normalization of the fft and ifft below (see costant c) is needed to make one operator 
%            pass the dot product test. In other words, the normalization is to guarantee that option = -1 is the
%            adjoint (transpose or conjugate operator) of option = 1. Please, see Claerbout's book Processing versus
%            inversion where he discusses the dot product test and the utilization of CG (or IRLS in my case) ''on the flight''.
%
% ------------------------------------------
% Mauricio D Sacchi, Natal, August 2015
% msacchi@ualberta.ca
% ------------------------------------------
%
% 


if  NDim == 2; 
 c = K(1);
  if (option == -1)
   aux = T.*input;
   output = fft(aux,K(1));
  end;
  if (option == 1); 
   aux = c*ifftn(input);
   output = T.*aux(1:N(1));
  end
end

if  NDim == 3; 
 c = K(1)*K(2);
  if (option == -1)
   aux = T.*input;
   output = fftn(aux,[K(1),K(2)]);
  end;
  if (option == 1); 
   aux = c*ifftn(input);
   output = T.*aux(1:N(1),1:N(2));
  end
end

if  NDim == 4; 
 c = K(1)*K(2)*K(3);
  if (option == -1)
   aux = T.*input;
   output = fftn(aux,[K(1),K(2),K(3)]);
  end;
  if (option == 1); 
   aux = c*ifftn(input);
   output = T.*aux(1:N(1),1:N(2),1:N(3));
  end
end

if  NDim == 5; 
 c = K(1)*K(2)*K(3)*K(4);
  if (option == -1)
   aux = T.*input;
   output = fftn(aux,[K(1),K(2),K(3),K(4)]);
  end;
  if (option == 1); 
   aux = c*ifftn(input);
   output = T.*aux(1:N(1),1:N(2),1:N(3),1:N(4));
  end
end

return
  
