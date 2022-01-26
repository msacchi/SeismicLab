function M=matrandorth(n, tol)
%MATRANDORTH Generates random n x n orthogonal real matrix.
%
%   M = MATRANDORTH(N) generates a random N x N orthogonal real matrix.
%
%   M = MATRANDORTH(M,TOL) explicitly specifies a threshold value, TOL,
%   that measures linear dependence of a newly formed column with the
%   existing columns. Defaults to 1e-6. 
%
%   In this version the generated matrix distribution *is* uniform over the
%   manifold O(n) w.r.t. the induced R^(n^2) Lebesgue measure, at a slight
%   computational overhead (randn + normalization, as opposed to rand ). 
% 
%   NOTE: This code is renamed from RandOrthMat by Ofek Shilon.
%   https://www.mathworks.com/matlabcentral/fileexchange/11783-randorthmat
%
%   (c) Ofek Shilon, 2006.
%
%   This code is *not* copyrighted by Sandia, but it is distributed with:
%
%   See also MATRANDNORM, MATRANDCONG, CREATE_PROBLEM, CREATE_GUESS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

    if exist('tol','var')
        if (tol >= 1)
            ncols = tol;
            tol = 1e-6;
        else
            ncols = n;
        end
    else
        tol=1e-6;
        ncols = n;
    end

    
    M = zeros(n); % prealloc
    
    % gram-schmidt on random column vectors
    
    vi = randn(n,1);  
    % the n-dimensional normal distribution has spherical symmetry, which implies
    % that after normalization the drawn vectors would be uniformly distributed on the
    % n-dimensional unit sphere.

    M(:,1) = vi ./ norm(vi);
    
    for i=2:n
	  nrm = 0;
	  while nrm<tol
		vi = randn(n,1);
		vi = vi -  M(:,1:i-1)  * ( M(:,1:i-1).' * vi )  ;
		nrm = norm(vi);
	  end
	  M(:,i) = vi ./ nrm;

    end %i
        
    M = M(:,1:ncols);
    
end  % RandOrthMat
    
