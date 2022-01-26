function out = rqrd(in,p)
%RQRD: ramdomized QR decomposition algorithm used to 
%      compute the low rank approximation of a matrix  
%
%  out = rqrd(in,p)
%
%  IN   in:      input matrix  
%       p:       desired rank 
%
%  OUT  out:     reduced-rank matrix
% 
%  Reference: Liberty et al., 2007, Randomized algorithms for the low-rank approximation 
%             of matrices, PNAS,104 (51) 20167-20172
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.

[~,ny]=size(in);

omega=rand(ny,p);

y=in*omega;

[U,~]=qr(y,0);

out=U*U'*in;

return


