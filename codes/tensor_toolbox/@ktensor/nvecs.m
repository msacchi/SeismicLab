function u = nvecs(X,n,r,opts)
%NVECS Compute the leading mode-n vectors for a ktensor.
%
%   U = NVECS(X,N,R) computes the R leading eigenvalues of Xn*Xn'
%   (where Xn is the mode-N matricization of X), which provides
%   information about the mode-N fibers. In two-dimensions, the R
%   leading mode-1 vectors are the same as the R left singular vectors
%   and the r leading mode-2 vectors are the same as the R right
%   singular vectors. By default, this method computes the top R
%   eigenvectors of the matrix Xn*Xn' using EIGS.
%
%   Examples
%   K = ktensor(@rand, [3 4 5], 2);
%   nvecs(K, 3, 1) %<--The largest eigenvector of the 3-mode matricization
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','nvecs_doc.html')))">Documentation page for n-vecs</a>
%
%   See also KTENSOR, TENMAT, EIGS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%   HIDDEN OPTIONS...
%
%   U = NVECS(X,n,r,OPTS) specifies additional options:
%   OPTS.eigsopts: options passed to the EIGS routine [struct('disp',0)]
%   OPTS.flipsign: make each column's largest element positive [true]
%

if ~exist('opts','var')
    opts = struct;
end

if isfield(opts,'eigsopts')
    eigsopts = opts.eigsopts;
else
    eigsopts.disp = 0;
end

% Compute Xn * Xn' excluding the nth factor
M = X.lambda * X.lambda';
for i = 1:ndims(X)
    if i == n, continue, end;
    M = M .* (X.u{i}' * X.u{i});
end

% Compute Xn * Xn'
Y = X.u{n} * M * X.u{n}';

[u,d] = eigs(Y, r, 'LM', eigsopts);

if isfield(opts,'flipsign') 
    flipsign = opts.flipsign;
else
    flipsign = true;
end
    
if flipsign
    % Make the largest magnitude element be positive
    [val,loc] = max(abs(u));
    for i = 1:r
        if u(loc(i),i) < 0
            u(:,i) = u(:,i) * -1;
        end
    end
end
