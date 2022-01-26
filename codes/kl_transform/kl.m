function [D_out,s] = kl(D_in,P);  
%KL: Karhunen-Loeve filtering of seismic data.
%    The Karhunen-Loeve transform can be used to enhance
%    lateral coherence of seismic events, this transform 
%    works quite well with NMO corrected CMP gathers, or
%    common offset sections.
%
%  [D_out,s] = kl(D_in,P)
%
%  IN   D_in:  data (the data matrix)
%       P:     number of eigenvectors to reconstruct the
%              data
%
%  OUT  D_out: output data after kl filtering
%       s:     the P largest eigen-values of the covariance
%              matrix in descending order
%
%  References: Jones, I.F. and Levy, S. 1987. Signal-to-noise ratio enhancement 
%              in multichannel seismic data via the  
%              Karhunen-Loeve transform, Geophysical Prospecting 35,12-32. 
%    
%              Al-Yahya Kamal, 1991, Application of the partial KL transform
%              to suppress random noise in seismic sections,
%              Geophysical Prospecting 39, 77-93.
%
%  Example:
%
%     d = flat_events(1); do = kl(d,3); wigb([d,do]);
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.





 [nt,nh]=size(D_in);
 
 D_out=zeros(nt,nh);

 R = D_in'*D_in;                   % Data covariance 

 OPTS.disp=0;

 [U,S]=eigs(R,P,'LM',OPTS);        % Eigen-decomposition of R

 U=U(:,1:P);             

 D_out=(U*U'*D_in')';

 s = diag(S);           
                        
 return


