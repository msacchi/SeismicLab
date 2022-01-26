function X = redistribute(X,mode)
%REDISTRIBUTE Distribute lambda values to a specified mode. 
%
%   K = REDISTRIBUTE(K,N) absorbs the weights from the lambda vector
%   into mode N. The lambda vector is then set to all ones.
%
%   Examples
%   K = ktensor([2; 4], ones(3,2), ones(4,2), ones(2,2));
%   redistribute(K, 3) %<--Weight vector is absorbed into factor matrix
%
%   See also KTENSOR, NORMALIZE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


for r = 1:length(X.lambda)
    X.u{mode}(:,r) = X.u{mode}(:,r) * X.lambda(r);
    X.lambda(r) = 1;
end
