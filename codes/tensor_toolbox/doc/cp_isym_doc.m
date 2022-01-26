%% Implicit Symmetric CP Decomposition for Symmetric K-Tensors
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="cp.html">CP Decompositions</a> 
% &#62;&#62; <a href="cp_isym_doc.html">CP-ISYM</a>
% </p>
% </html>
%
% The function |cp_isym| computes the symmetric CP
% decomposition of a symmetric tensor that is in symmetric k-tensor format,
% meaning that it is the sum of symmetric outer products and stored as a
% |symktensor| object.
% The decomposition is described in the following reference:
%
% * S. Sherman, T. G. Kolda, Estimating Higher-Order Moments Using
%   Symmetric Tensor Decomposition, SIAM J. Matrix Analysis and 
%   Applications, 41:1369-1387, 2020,
%   <https://doi.org/10.1137/19m1299633>

%% Requirements
% This code requires an optimization solver. Our examples use the L-BFGS-B
% optimization method; see <opt_options_doc.html Optimization Options> for
% details of installation and usage as well as other options for
% optimization solvers.


%% Create a sample problem based on a Gaussian mixture model. 
% Here we used the third-order moment (d=3), for a problem with random
% variables of dimension n=100. We take a mixture of r=5 Gaussians and
% collect p=750 observations.
d = 3; 
n = 100;
r = 5;
p = 250*r;

%%
% We construct the means so that they are collinear, i.e., the cosine of
% the angle between any two vectors is 0.5. This is a more difficult
% problem than purely random vectors which are close to orthogonal.
rng('default');
collinearity = 0.5;
Atrue = matrandcong(n,r,collinearity);

%%
% Each Gaussian gets equal weight in the mixture.
wghts = 1/r * ones(r,1);

%%
% The `true` solution is based on the means and their prevalence in the
% mixture.
Mtrue = symktensor(wghts,Atrue,d);
Mtrue = normalize(Mtrue);

%%
% The data tensor is based on noisy observations from the mixture model.
% The |idx| determines which Gaussian is sampled. And the samples have
% noise as specified by |stddev|. The result is assembled into an
% observation tensor that in |symktensor| format.
idx = randsample(r,p,true,wghts);
stddev = 1e-3;
noise = stddev * randn(n,p);
V = Atrue(:,idx) + noise;
X = symktensor(1/p*ones(p,1),V,d);
X = normalize(X,0);

%% Recommend multiple runs of CP-ISYM 
% Since this is a non-convex optimization problem, it's critically
% important to do multiple runs of the optimzation problem. In the code
% below, we run 10 times. Not every run will converge to the global minimum. 
% We take solution with the lowest function value, which should ideally
% yield the best match to the true solution. For our artificial problem, we
% can actually verify the score of how well the computed solution matches
% with the ideal solution.
rng(7)
nruns = 10;
M = cell(nruns,1);
info = cell(nruns,1);
for i = 1:nruns
    fprintf('------------- Run %d of %d --------------\n', i, nruns);
    [M{i},info{i}] = cp_isym(X,r,'Xnormsqr','exact','printitn',0);
    fprintf('Final F: %.2e, Score: %.2f\n', info{i}.f, score(M{i},Mtrue));
    if i == 1
        ifbest = i;
        fbest = info{i}.f;
        sfbest = score(M{i},Mtrue);
    elseif info{i}.f < fbest
        ifbest = i;
        fbest = info{i}.f;
        sfbest = score(M{i},Mtrue);
    end
end

%%
fprintf('------------- Best of %d runs --------------\n', nruns);
fprintf('Run %d -> Final F: %.2e, Score: %.2f\n', ifbest, fbest, sfbest);

%% Options for CP-ISYM with L-BFGS-B (default) optimization solver
% The CP solution should produce something close to |Mtrue|, but it will
% not be exact because of the randomness in building the observation
% tensor. Here we show what happens when the true solution is provided as
% an initial guess. The function should converge to zero, and the
% similarlity score should be close to one. The method defaults to
% |'lfbgsb'|.

[M,info] = cp_isym(X,r,'Xnormsqr','exact','init',Mtrue);
fprintf('\n\t\t*** Final F: %.2e, Score: %.2f ***\n\n', info.f, score(M,Mtrue));

%% 
% It is also possible to specify the initial guess as factor matrix 
[M,info] = cp_isym(X,r,'Xnormsqr','exact','init',Mtrue);
fprintf('\n\t\t*** Final F: %.2e, Score: %.2f ***\n\n', info.f, score(M,Mtrue));

%%
% So far we've been using the exact objective function which is close to
% zero for the optimal zero. We can, however, ignore the constant ||X||^2
% term in objective function, which saves some preprocessing cost. The main
% disadvantage of ignoring the constant term is that the final function
% value will no longer be near zero. 
[M,info] = cp_isym(X,r,'Xnormsqr',0,'init',Atrue,'printitn',1);
fprintf('\n\t\t*** Final F (adjusted): %.2e, Score: %.2f ***\n\n', ...
    info.f + norm(X)^2, score(M,Mtrue));

%% 
% By default, the method uses the randomized range finder (rrf) for the
% initial guess. This is the default initial guess and works well.
rng(1)
[M,info] = cp_isym(X,r,'Xnormsqr','exact','printitn',5,'init','rrf');
fprintf('\n\t\t*** Final F: %.2e, Score: %.2f ***\n\n', info.f, score(M,Mtrue));


%% 
% If we run instead with a random Gaussian initial guess, we are unlikey to
% converge. 
rng(77)
[M,info,M0] = cp_isym(X,r,'Xnormsqr','exact','printitn',5,'init',@randn);
fprintf('\n\t\t*** Final F: %.2e, Score: %.2f ***\n\n', info.f, score(M,Mtrue));


%%
% Sometimes a random initialization does a little better if the convergence
% tolerances are tightened and the maximum iterations is increased. We also
% increase the memory in Limited-memory BFGS to see if that helps. It gets
% a little better, but still no where near the solution with the randomized
% range finder. 
[M,info] = cp_isym(X,r,'Xnormsqr','exact','printitn',5,'init',M0,...
    'ftol',1e-14,'gtol',1e-8,'maxiters',10000,'m',25);
fprintf('\n\t\t*** Final F: %.2e, Score: %.2f ***\n\n', info.f, score(M,Mtrue));

%% Options for CP-ISYM with L-BFGS optimization solver from Poblano
rng(1)
[M,info] = cp_isym(X,r,'Xnormsqr','exact','printitn',5,'method','lbfgs');
fprintf('\n\t\t*** Final F: %.2e, Score: %.2f ***\n\n', info.f, score(M,Mtrue));

%% Options for CP-ISYM with Quasi-Newton optimization solver from Optimization Toolbox
rng(1)
[M,info] = cp_isym(X,r,'Xnormsqr','exact','printitn',5,'method','fminunc');
fprintf('\n\t\t*** Final F: %.2e, Score: %.2f ***\n\n', info.f, score(M,Mtrue));

%% Options for CP-ISYM with Adam stochastic optimization solver
% When the number of observations is very large, a stochastic method may be
% faster. Our example is very small (p=1250), but we use it nonetheless
% just to demonstrate the capabilities of the stochastic method. The method
% is called by specifying the |'method'| to be |'adam'|. 

rng(5)
[M,info] = cp_isym(X,r,'method','adam');
fprintf('\n\t\t*** Final F: %.2e, Score: %.2f ***\n\n', f_implicit(M,X,norm(X)^2), score(M,Mtrue));

%%
% Since this problem is small, we can compute the function value exactly.
rng(5)
[M,info] = cp_isym(X,r,'method','adam','fsamp','exact','Xnormsqr','exact');
fprintf('\n\t\t*** Final F: %.2e, Score: %.2f ***\n\n', f_implicit(M,X,norm(X)^2), score(M,Mtrue));

%%
% We can also change the number of samples in each stochastic gradient
% (|'gsamp'|) and correspondingly increase the number of iterations per
% epoch (|'subiters'|).
rng(5)
[M,info] = cp_isym(X,r,'method','adam','fsamp','exact','Xnormsqr','exact',...
    'gsamp',10,'subiters',250);
fprintf('\n\t\t*** Final F: %.2e, Score: %.2f ***\n\n', f_implicit(M,X,norm(X)^2), score(M,Mtrue));

%%
% We can further improve by increasing the number of times we decrease the
% step length (|'maxfails'|). 
rng(5)
[M,info] = cp_isym(X,r,'method','adam','fsamp','exact','Xnormsqr','exact',...
    'gsamp',10,'subiters',250,'maxfails',4);
fprintf('\n\t\t*** Final F: %.2e, Score: %.2f ***\n\n', f_implicit(M,X,norm(X)^2), score(M,Mtrue));

%%
% We also change the initial learning rate via 'rate'.
rng(5)
[M,info] = cp_isym(X,r,'method','adam','fsamp','exact','Xnormsqr','exact','rate',1e-3);
fprintf('\n\t\t*** Final F: %.2e, Score: %.2f ***\n\n', f_implicit(M,X,norm(X)^2), score(M,Mtrue));
