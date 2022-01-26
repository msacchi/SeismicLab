function A_out = rand_svd(A,k,l,i);
%RAND_SVD: randomized SVD to compute the low-rank approximation
%          of a matrix 
%
%  A_out = rand_svd(A,K)
%
%  IN   A:      input matrix 
%       k:      desired rank 
%       l = 2*k seems to work
%       i = 2   also is a good option
%
%  OUT  A_out:  reduced-rank matrix 
% 
%  Reference: Liberty et al., 2007, Randomized algorithms for the low-rank approximation 
%             of matrices, PNAS,104 (51) 20167-20172
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.
%

[m,n] = size(A);

GG = rand([l,m]); % Random l x m R matrix.

AA = 1;

for ii=1:i
    AA = AA*(A*A');
end

RR = GG*AA*A;

[UUq,EEq,VVq] = svd(RR',0);

QQ = UUq(:,1:k);

TT = A*QQ;

[UU,EE,W] = svd(TT,0);

VV = QQ*W;

A_out = UU*EE*VV';

return;
