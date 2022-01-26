function [F,G] = fg(Model,Data)
%FG Master objective function for optimization of symmetric Kruskal model.
%
%   [F,G] = FG(MODEL,DATA) computes the function (F) and gradient (G) for
%   the given MODEL and DATA.  The input DATA should be computed by calling
%   the FG_SETUP function. The function value F is a scalar and the
%   gradient G is returned as a column vector of length Q = P*(N+1) (or
%   Q=P*N if 'nolambda' is true in FG_SETUP). Here, P is the rank of the
%   decomposition and N is the number of modes in MODEL.
%
%   Example:
%   A = symmetrize(tenrand(2,2,2)); 
%   P = 3; 
%   model = symktensor(P,A); %<- Create random symtensor of size 2x2x2
%   data = fg_setup(model,A,'unique',false);
%   [f,g] = fg(model,data);
%   f - norm(A - full(model)).^2 % Should be zero or close to it
%
%   See also symktensor, cp_sym, symktensor/fg_setup, symktensor/tovec.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%% F/G computation

% Extract model
X = Model.u;
lambda = Model.lambda;

% Extract options
fgopts = Data.fgopts;

if fgopts.fast

    % Extract data
    M = Data.M;
    N = Data.N;
    P = Data.P;
    A = Data.A;
    
    % Precomputation
    XtX = X'*X;
    XtXMM = XtX.^(M-1);
    XtXM = XtXMM .* XtX;
    %XtXM = XtX.^M;
    Y = zeros(size(X));
    Z = zeros(P,1);
    for p = 1:P
        Y(:,p) = ttsv(A,X(:,p),-1); % Ax^(m-1)
        Z(p) = dot(X(:,p),Y(:,p)); % Ax^m
    end
    
    % Compute F
    F1 = Data.normAsqr;
    F2 = -2 * dot(lambda,Z);
    F3 = lambda' * XtXM * lambda;
    F = F1 + F2 + F3;
    
    % Compute G
    
    % Compute G wrt Lambda
    if ~fgopts.nolambda
        G.lambda = -2 * Z + 2 * XtXM * lambda;
    end
    
    % Compute G wrt X
    G.X = -2 * M * Y * diag(lambda) + 2 * M * X*diag(lambda)*XtXMM*diag(lambda);
        
    % Convert G to vector
    G = stackgrad(G,fgopts.nolambda);
    
    
else

    % Extract data
    M = Data.M;
    N = Data.N;
    P = Data.P;
    Q = Data.Q;
    %R = Data.R;
    avals = Data.avals;
    I = Data.I;
    C = Data.C;
    W = Data.W;
    onesidx = Data.onesidx;
    zerosidx = Data.zerosidx;
    %lb = Data.lb;
    %ub = Data.ub;
    
    % --- Compute F/G
    
    % Compute diffs
    foo = X(I,:);
    foo = reshape(foo, [Q M P]);
    xprods = squeeze(prod(foo,2)); % q x p matrix
    vals = xprods * lambda;
    diffs = avals - vals;
    
    % Compute F
    F = sum(W.*(diffs.^2));
    
    % Compute G wrt Lambda
    if ~fgopts.nolambda
        for k = 1:P
            G.lambda(k,1) = -2 * sum( W .* diffs .* xprods(:,k) );
        end
    end
    
    % Compute G wrt X
    for n = 1:N
        bar = foo;
        bar = reshape(bar,Q*M,P);
        bar(onesidx{n},:) = 1;
        bar = reshape(bar,Q,M,P);
        bar(zerosidx{n},:,:) = 0;
        xprodsjm1 = squeeze(prod(bar,2));
        for k = 1:P
            G.X(n,k) = -2 * lambda(k) * sum( C(:,n) .* W .* diffs .* xprodsjm1(:,k) );
        end
    end
    
    % Convert G to vector
    G = stackgrad(G,fgopts.nolambda);
end
% --- Penalties ---

% Norm weight
if (fgopts.l2weight > 0)
    tmp = bsxfun(@dot,X,X)-1;
    ftmp = fgopts.l2weight * sum(tmp.^2);
    gtmp.lambda = zeros(P,1);
    gtmp.X = fgopts.l2weight * 4 * bsxfun(@mtimes,X,tmp);
    F = F + ftmp;
    G = G + stackgrad(gtmp);
end

% Lambda weight
if (fgopts.l1weight > 0)
    alpha = fgopts.l1param;
    lambda_alpha = (1/alpha) * (log(1 + exp(-alpha*lambda)) + log(1+exp(alpha*lambda)));
    ftmp = fgopts.l1weight * sum(lambda_alpha);
    grad_alpha = 1./(1 + exp(-alpha*lambda)) - 1./(1 + exp(alpha*lambda));
    gtmp.lambda = fgopts.l1weight * grad_alpha;
    gtmp.X = zeros(N,P);
    F = F + ftmp;
    G = G + stackgrad(gtmp);
end


function g = stackgrad(G,nolambda)
if exist('nolambda','var') && nolambda
    g = reshape(G.X,[],1);
else    
    g = [G.lambda; reshape(G.X,[],1)];
end
