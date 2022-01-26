function n = ncomponents(t)
%NCOMPONENTS Number of components for a ktensor.
%
%   NCOMPONENTS(T) returns the number of components in the ktensor T.
%
%   X = ktensor(ones(4,1), rand(2,4), randn(3,4), randi(5,4,4));
%   ncomponents(X) %<--Returns 4
%
%   See also KTENSOR
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



n = length(t.lambda);
