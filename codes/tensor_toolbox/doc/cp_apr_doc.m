%% Alternating Poisson Regression for fitting CP to sparse count data
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="cp.html">CP Decompositions</a> 
% &#62;&#62; <a href="cp_apr_doc.html">CP-APR</a>
% </p>
% </html>
%
% Reference: E. C. Chi, T. G. Kolda, On Tensors, Sparsity, and Nonnegative Factorizations,
% SIAM J. Matrix Analysis and Applications, 33:1272-1299, 2012, https://doi.org/10.1137/110859063.
%

%% Set up a sample problem
% We follow the general procedure outlined by Chi and Kolda (2013).

% Pick the size and rank
sz = [100 80 60];
R = 5;

% Generate factor matrices with a few large entries in each column; this
% will be the basis of our soln.
A = cell(3,1);
for n = 1:length(sz)
    A{n} = rand(sz(n), R);
    for r = 1:R
        p = randperm(sz(n));
        nbig = round( (1/R)*sz(n) );
        A{n}(p(1:nbig),r) = 100 * A{n}(p(1:nbig),r);
    end
end
lambda = rand(R,1);
S = ktensor(lambda, A);
S = normalize(S,'sort',1);

% Create sparse test problem based on provided solution. 
nz = prod(sz) * .05;
info = create_problem('Soln', S, 'Sparse_Generation', nz);

% Extract data and solution
X = info.Data;
M_true = info.Soln;

%% Call CP-APR

% Compute a solution
M = cp_apr(X, R, 'printitn', 10);

% Score the solution
factor_match_score = score(M, M_true, 'greedy', true)
