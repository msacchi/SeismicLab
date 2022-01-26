%% Alternating randomized least squares for CP Decomposition
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="cp.html">CP Decompositions</a> 
% &#62;&#62; <a href="cp_arls_doc.html">CP-ARLS</a>
% </p>
% </html>
%
% The function |cp_arls| computes an estimate of the best rank-R CP model
% of a tensor X using alternating _randomized_ least-squares algorithm.
% The input X must be a (dense) |tensor|. The output CP model is a
% |ktensor|.  The CP-ARLS method is described in the following reference:
%
% * C. Battaglino, G. Ballard, T. G. Kolda.
%   A Practical Randomized CP Tensor Decomposition.
%   SIAM Journal on Matrix Analysis and Applications 39(2):876-901, 2018.
%   <https://doi.org/10.1137/17M1112303>
%% Set up a sample problem
% We set up an especially difficult and somewhat large sample problem that
% has high collinearity (0.9) and 1% noise. This is an example where the
% randomized method will generally outperform the standard method.
sz = [200 300 400];
R = 5;
ns = 0.01;
coll = 0.9;

info = create_problem('Size', sz, 'Num_Factors', R, 'Noise', ns, ...
    'Factor_Generator', @(m,n) matrandcong(m,n,coll), ...
    'Lambda_Generator', @ones);

% Extract data and solution
X = info.Data;
M_true = info.Soln;

%% Running the CP-ARLS method
% Running the method is essentially the same as using CP-ALS, feed the data
% matrix and the desired rank. Note that the iteration is of the form NxN
% which is the number of epochs x the number of iterations per epoch. The
% default number of iterations per epoch is 50. At the end of each epoch,
% we check the convergence criteria. Because this is a randomized method,
% we do not achieve strict decrease in the objective function. Instead, we
% look at the number of epochs without improvement (newi) and exit when
% this crosses the predefined tolerance (|newitol|), which defaults to 5.
% It is important to note that the fit values that are reported are
% approximate, so this is why it is denoted by |f~| rather than just |f|.

tic
[M1, ~, out1] = cp_arls(X,R);
time1 = toc;
scr1 = score(M1,M_true);
fprintf('\n*** Results for CP-ARLS (with mixing) ***\n');
fprintf('Time (secs): %.3f\n', time1)
fprintf('Score (max=1): %.3f\n', scr1);

%% Speed things up by skipping the initial mixing
% The default behavior is to mix the data in each mode using an FFT and
% diagonal random +/-1 matrix. This may add substantial preprocessing time,
% though it helps to ensure that the method converges. Oftentimes, such as
% with randomly-generated data, the mixing is not necessary. 

tic
[M2, ~, out2] = cp_arls(X,R,'mix',false);
time2 = toc;
scr2 = score(M2,M_true);

fprintf('\n*** Results for CP-ARLS (no mix) ***\n');
fprintf('Time (secs): %.3f\n', time2)
fprintf('Score (max=1): %.3f\n', scr2);

%% Comparing with CP-ALS
% CP-ALS may be somewhat faster, especially since this is a relatively
% small problem, but it usually will not achieve as good of an answer in
% terms of the score.

tic; 
[M3, ~, out3] = cp_als(X,R,'maxiters',500,'printitn',10); 
time3 = toc;
scr3 = score(M3,M_true);
fprintf('\n*** Results for CP-ALS ***\n');
fprintf('Time (secs): %.3f\n', time3)
fprintf('Score (max=1): %.3f\n', scr3);

%% How well does the approximate fit do?
% It is possible to check the accuracy of the fit computation by having the
% code compute the true fit and the final solution, enabled by the
% |truefit| option.
[M4,~,out4] = cp_arls(X,R,'truefit',true);

%% Varying epoch size
% It is possible to vary that number of iterations per epoch. Fewer
% iterations means that more time is spent checking for convergence and it
% may also be harder to detect as an single iteration can have some
% fluctuation and we are actually looking for the overall trend. In
% contrast, too many iterations means that the method won't realize when it
% has converged and may spend too much time computing.

%%
tic
M = cp_arls(X,R,'epoch',1,'newitol',20);
toc
fprintf('Score: %.4f\n',score(M,M_true));

%%
tic
M = cp_arls(X,R,'epoch',200,'newitol',3,'printitn',2);
toc
fprintf('Score: %.4f\n',score(M,M_true));

%% Set up another sample problem
% We set up another problem with 10% noise, but no collinearity.
sz = [200 300 400];
R = 5;
ns = 0.10;

info = create_problem('Size', sz, 'Num_Factors', R, 'Noise', ns, ...
    'Factor_Generator', @rand,'Lambda_Generator', @ones);

% Extract data and solution
X = info.Data;
M_true = info.Soln;

%% Terminating once a desired fit is achieved
% If we know the noise level is 10%, we would expect a fit of 0.90 at best.
% So, we can set a threshold that is close to that and terminate as soon as
% we achieve that accuracy. Since detecting convergence is hard for a
% randomized method, this can lead to speed ups. However, if the fit is not
% high enough, the accuracy may suffer consequently.
M = cp_arls(X,R,'newitol',20,'fitthresh',0.895,'truefit',true);
fprintf('Score: %.4f\n',score(M,M_true));

%% Changing the number of function evaluation samples
% The function evaluation is approximate and based on sampling the number
% of entries specified by |nsampfit|. If this is too small, the samples
% will not be accurate enough. If this is too large, the computation will
% take too long. The default is $2^14$, which should generally be
% sufficient.  It may sometimes be possible to use smaller values. The same
% sampled entries are used for every convergence check --- we do not
% resample to check other entries.
M = cp_arls(X,R,'truefit',true,'nsampfit',100);
fprintf('Score: %.4f\n',score(M,M_true));

%% Change the number of sampled rows in least squares solve
% The default number of sampled rows for the least squares solves is
% |ceil(10*R*log2(R))|. This seemed to work well in most tests, but this can
% be varied higher or lower. For R=5, this means we sample 117 rows per
% solve. The rows are different for every least squares problem. Let's see
% what happens if we reduce this to 10.

M = cp_arls(X,R,'truefit',true,'nsamplsq',10);
fprintf('Score: %.4f\n',score(M,M_true));

%%
% What if we use 25?
M = cp_arls(X,R,'truefit',true,'nsamplsq',25);
fprintf('Score: %.4f\n',score(M,M_true));


