function nrm = norm(X)
%NORM Norm of a ttensor.
%
%   See also TTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if prod(size(X)) > prod(size(X.core))
    V = cell(ndims(X),1);
    for n = 1:ndims(X)
        V{n} = X.u{n}'*X.u{n};
    end
    Y = ttm(X.core,V);
    tmp = innerprod(Y, X.core);
    nrm = sqrt(tmp);
else
    nrm = norm(full(X));
end
