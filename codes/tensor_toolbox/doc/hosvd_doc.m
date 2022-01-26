%% Computing Tucker via the HOSVD
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="tucker.html">Tucker Decompositions</a> 
% &#62;&#62; <a href="hosvd_doc.html">HOSVD</a>
% </p>
% </html>
%

%% Higher-order Singular Value Decomposition (HOSVD) and Sequentially-truncased HOSVD (ST-HOSVD)
% The HOSVD computes a Tucker decomposition of a tensor via a simple
% process. For each mode k, it computes the r_k leading left singular
% values of the matrix unfolding and stores those as factor matrix U_k. 
% Then it computes a |ttm| of the original tensor and all the factor matrices to
% yield the core of size r_1 x r_2 x ... x r_d. The core and factor
% matrices are used to form the |ttensor|. 
% The values of r_k that lead to a good approximation can be computed
% automatically to yield a specified error tolerance; this is recommended
% and the default in our code.
% The ST-HOSVD is an improvement on the HOSVD that does a TTM in _each_ mode
% before moving on to the next mode. This has the advantage of shrinking
% the tensor at each step and reducing subsequent computations. ST-HOSVD is the
% default in the |hosvd| code.
%
%
% * L. R. Tucker, Some mathematical notes on three-mode factor analysis,
%   Psychometrika, 31:279-311, 1966, http://dx.doi.org/10.1007/BF02289464
% * L. D. Lathauwer, B. D. Moor and J. Vandewalle, A multilinear singular
%   value decomposition, SIAM J. Matrix Analysis and Applications,
%   21(4):1253-1278, 2000, http://dx.doi.org/10.1137/S0895479896305696  
% * N. Vannieuwenhoven, R. Vandebril and K. Meerbergen, A New Truncation
%   Strategy for the Higher-Order Singular Value Decomposition, SIAM J.
%   Scientific Computing, 34(2):A1027-A1052, 2012,
%   http://dx.doi.org/10.1137/110836067    
%

%% Simple example of usage

% Create random 50 x 40 x 30 tensor with 5 x 4 x 3 core
info = create_problem('Type','Tucker','Num_Factors',[5 4 3],'Size',[50 40 30],'Noise',0.01);
X = info.Data;

% Compute HOSVD with desired relative error = 0.1
T = hosvd(X,0.1);

% Check size of core
coresize = size(T.core)

% Check relative error
relerr = norm(X-full(T))/norm(X)

%% Generate a core with different accuracies for different sizes
% We will create a core tensor that has is nearly block diagonal. The
% blocks are expontentially decreasing in norm, with the idea that we can
% pick off one block at a time as we increate the prescribed accuracy of
% the HOSVD. To do this, we use |tenrandblk|.

% Block sizes (need not be cubic). Number of rows is the number
% of levels and number of columns is the order of the tensor.
bsz = [3 2 1; 2 2 2; 2 3 4];

% Squared norm of each block. Must be length L and sum to <= 1
bsn = [.9 .09 .009]';

% Create core tensor with given block structure and norm 1
G = tenrandblk(bsz,bsn,true);

%%
fprintf('Size of G: %s\n', tt_size2str(size(G)));

%% Generate data tensor with core described above
% We take the core G and embed into into a larger tensor X by using
% orthogonal transformations. The true rank of this tensor is equal to the
% size of G.

% Size of X
xsz = [20 20 20];

% Create orthogonal matrices
U = cell(3,1);
for k = 1:3
    V = matrandorth(xsz(k));
    U{k} = V(:,1:size(G,k));
end

% Create X
X = full(ttensor(G,U));

% The norm should be unchanged
fprintf('||X||=%f\n',norm(X));

%% Compute (full) HOSVD 
% We compute the ST-HOSVD using the |hosvd| method. We specify the
% tolerance to close to machine precision. Ideally, it finds a core that is
% the same size as G. 

fprintf('ST-HOSVD...\n');
T = hosvd(X,2*sqrt(eps));

%% Compute low-rank HOSVD approximation
% The norm squared of the first two blocks of G is 0.99, so specifying an
% error of 1e-2 should yield a core of size 4 x 4 x 3.  However, the
% conservative nature of the algorithm means that it may pick something
% larger. We can compensate by specifying a larger tolerance.

% Using 1e-2 exactly is potentially too conservative...
fprintf('Result with tol = sqrt(1e-2):\n');
T = hosvd(X, sqrt(1e-2),'verbosity',1);

% But a small multiple (i.e., |ndims(X)|) usually works...
fprintf('Result with tol = sqrt(3e-2):\n');
T = hosvd(X, sqrt(3e-2),'verbosity',1);

%%
% Similarly, lhe norm squared of the first block of G is 0.9, so specifying
% an error of 1e-1 should result in a core of size 3 x 2 x 1.  

% Using 1e-1 exactly is potentially too conservative...
fprintf('Result with tol = sqrt(1e-1):\n');
T = hosvd(X, sqrt(1e-1),'verbosity',1);

% But a small multiple (i.e., |ndims(X)|) usually works...
fprintf('Result with tol = sqrt(3e-1):\n');
T = hosvd(X, sqrt(3e-1),'verbosity',1);

%% Verbosity - Getting more or less information.
% Setting the verbosity to zero suppresses all output.
% Cranking up the verbosity gives some insight into the decision-making
% process...

% Example 1
T = hosvd(X, sqrt(3e-1),'verbosity',10);

%%
% Example 2
T = hosvd(X, sqrt(3*eps),'verbosity',10);

%% Specify the ranks
% If you know the rank  you want, you can specify it. But there's no
% guarantee that it will satisfy the specified tolerance. In such cases,
% the method will throw a warning.

% Rank is okay
T = hosvd(X,sqrt(3e-1),'ranks',bsz(1,:));

% Rank is too small for the specified error
T = hosvd(X,sqrt(3e-1),'ranks',[1 1 1]);

% But you can set the error to the tensor norm to make the warning go away
T = hosvd(X,norm(X),'ranks',[1 1 1]);

%% Specify the mode order
% It's also possible to specify the order of the modes. The default is
% 1:ndims(X).
T = hosvd(X,sqrt(3e-1),'dimorder',ndims(X):-1:1);

%% Generate bigger data tensor with core described above
% Uses the same procedure as before, but now the size is bigger.

% Size of Y
ysz = [100 100 100];

% Create orthogonal matrices
U = cell(3,1);
for k = 1:3
    V = matrandorth(ysz(k));
    U{k} = V(:,1:size(G,k));
end

% Create Y
Y = full(ttensor(G,U));

%% ST-HOSVD compared to HOSVD
% The answers are essentially the same for the sequentially-truncated HOSVD
% and the HOSVD...

fprintf('ST-HOSVD...\n');
T = hosvd(Y,2*sqrt(eps));
fprintf('HOSVD...\n');
T = hosvd(Y,2*sqrt(eps),'sequential',false);

%%
% But ST-HOSVD may be slightly faster than HOSVD for larger tensors.

fprintf('Time for 10 runs of ST-HOSVD:\n');
tic, for i =1:10, T = hosvd(Y,2*sqrt(eps),'verbosity',0); end; toc

fprintf('Time for 10 runs of HOSVD:\n');
tic, for i =1:10, T = hosvd(Y,2*sqrt(eps),'verbosity',0,'sequential',false); end; toc


