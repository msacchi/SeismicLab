function Y = symmetrize(X)
%SYMMETRIZE Symmetrize a ktensor X in all modes.
%
%   Y = symmetrize(X) will symmetrize a ktensor X with respect to all
%   modes so that Y is symmetric with respect to any permutation of
%   indices.
%
%   See also ISSYMMETRIC, SYMKTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%T. Kolda, June 2014

n = ndims(X);
sz = size(X);

% Check tensor dimensions for compatibility with symmetrization
if any(sz(2:end) ~= sz(1))
    error('Tensor is not cubic -- cannot be symmetrized');
end

% Distribute lambda evenly into factors
X = normalize(X,0);

lambda = X.lambda;
U = X.u;
U1 = U{1};

V = U1;
for i = 2:n

    Ui = U{i};
    
    for j = 1:size(U1,2);
        if dot( U1(:,j), Ui(:,j) ) < 0
            Ui(:,j) = -Ui(:,j);
            lambda(j) = -lambda(j);
        end
    end
    
    V = V + Ui;
end

V = V./ n;

% Odd-ordered tensors should not have any negative lambda values
if mod(ndims(X),2) == 1
    for j = 1:length(lambda)
        if lambda(j) < 0
            lambda(j) = -lambda(j);
            V(:,j) = -V(:,j);
        end
    end
end

Y = cell(n,1);
for i = 1:n
    Y{i} = V;
end
Y = ktensor(lambda,Y);

%Y = arrange(Y);




    



