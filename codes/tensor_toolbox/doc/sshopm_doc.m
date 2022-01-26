%% Shifted Symmetric Higher-order Power Method
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="eigen.html">Tensor Eigenproblem</a> 
% &#62;&#62; <a href="sshopm_doc.html">SSHOPM</a>
% </p>
% </html>
%
% The methods are described in the following paper:
%
%   * T. G. Kolda, J. R. Mayo, Shifted Power Method for Computing Tensor
%     Eigenpairs, SIAM J. Matrix Analysis and Applications, 32:1095-1124,
%     2011, http://dx.doi/org/10.1137/100801482 
%   * T. G. Kolda, J. R. Mayo, An Adaptive Shifted Power Method for
%     Computing Generalized Tensor Eigenpairs, SIAM J. Matrix Analysis and
%     Applications, 35:1563-1582, 2014, http://dx.doi.org/0.1137/140951758    
%
%
% Note that there is also a method for finding complex eigenpairs: |eig_sshopmc|.

%% Data tensor 
% From Example 1 in E. Kofidis and P. A. Regalia, On the best rank-1
% approximation of higher-order supersymmetric tensors, SIAM J. Matrix
% Anal. Appl., 23:863-884, 2002, http://dx.doi.org/10.1137/S0895479801387413.
A = tenzeros([3 3 3 3]);
A(perms([1 1 1 1])) = 0.2883;
A(perms([1 1 1 2])) = -0.0031;
A(perms([1 1 1 3])) = 0.1973;
A(perms([1 1 2 2])) = -0.2485;
A(perms([1 1 2 3])) = -0.2939;
A(perms([1 1 3 3])) = 0.3847;
A(perms([1 2 2 2])) = 0.2972;
A(perms([1 2 2 3])) = 0.1862;
A(perms([1 2 3 3])) = 0.0919;
A(perms([1 3 3 3])) = -0.3619;
A(perms([2 2 2 2])) = 0.1241;
A(perms([2 2 2 3])) = -0.3420;
A(perms([2 2 3 3])) = 0.2127;
A(perms([2 3 3 3])) = 0.2727;
A(perms([3 3 3 3])) = -0.3054;


%% 
% Check symmetry of result
issymmetric(A)

%% Call |eig_sshopm| with no shift (fails to converge)
rng('default')
[lambda, x, flag, it] = eig_sshopm(A, 'MaxIts', 100, 'Shift', 0,'Display',1);

%% Call |eig_sshopm| with automatic shift
rng('default')
[lambda, x, flag, it] = eig_sshopm(A, 'MaxIts', 100,'Display',1);

%% Call |eig_sshopm| with fixed shift 
rng('default')
[lambda, x, flag, it] = eig_sshopm(A, 'MaxIts', 100, 'Shift', 1,'Display',1);

%% Convert to |symtensor| object 
% Note that the |eig_sshopm| method actually expects a standard |tensor|
% object, but we just display the symmetric tensor here.
Asym = symtensor(A)
