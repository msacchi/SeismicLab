function C = mldivide(A,B)
%MLDIVIDE Slash left division for sparse tensors.
%
%   MlDIVIDE(A,B) is called for the syntax 'A \ B' when A is a scalar and B
%   is a sparse tensor. 
%
%   Example
%   X = sptenrand([4 3 2],5);
%   3 \ X
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isscalar(A)
    newsubs = B.subs;
    newvals = B.vals / A;
    if A == 0
        nansubs = setdiff(allsubs(A),newsubs,'rows');
        newsubs = [newsubs; nansubs];
        newvals = [newvals; repmat(NaN,size(nansubs,1),1)];
    end
    C = sptensor(newsubs,newvals,B.size);
    return;
end

error('MLDIVIDE only supports the scalar case for sparse tensors');
