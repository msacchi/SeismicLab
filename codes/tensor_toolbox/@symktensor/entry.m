function value = entry(model,I)
%ENTRY Extract a single entry from a symktensor.
%
%   V = ENTRY(M,I) returns the value of M at entry I. M is a symktensor and
%   I is a matrix with rows which are indexes to query.
%   This value is not stored explicitly and so must be computed.
%
%   Examples
%   S = symtenrand(3,4); <-- Declare a random symtensor of size [4,4,4]
%   M = cp_sym(S,2); <-- Decompose S into a symktensor with rank 2
%   I = [1,2,4; 4,4,4];  <-- Matrix of indices to query
%   entry(S,I) <-- Query elements [1,2,4] and [4,4,4]
%
%   See also SYMKTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%% Extract
lambda = model.lambda;
X = model.u;

%% Case I: I = row
if isrow(I)        
    value = dot(lambda,prod(X(I,:)));
    return;    
end


%% Case II: list of indices
q = size(I,1);
m = size(I,2);
p = size(X,2);
foo = X(I,:);
foo = reshape(foo, [q m p]);
% squeeze(foo(q,:,:)) =  X(I(q,:),:)
bar = squeeze(prod(foo,2));
value = bar * lambda;
