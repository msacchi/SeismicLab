function M = update(M,v)
%UPDATE Convert vector to symktensor.
%
%   M = UPDATE(M,V) overwrites some or all of the information in the
%   symktensor M with data from the vector V. Assume M is a symktensor of
%   dimension N with R components. If V is a vector of length (N+1)*R, then
%   the first R entries overwrite the weight vector M.lambda and the
%   remainder overwrites the factor matrix M.U. If V is a vector of length
%   N*R, then that data overwrites the factor matrix M.U and the M.lambda
%   is unchanged. This function is typically used in the context of 
%   optimization methods.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

[n,r] = size(M.u);

if size(v,1) == (n+1)*r
    M.lambda = v(1:r);
    M.u = reshape(v(r+1:end),[n,r]);
elseif size(v,1) == n*r
    M.u = reshape(v,[n r]);
else
    error('v is not the right size')
end
