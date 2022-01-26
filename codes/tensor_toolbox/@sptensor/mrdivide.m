function C = mrdivide(A,B)
%MRDIVIDE Slash right division for sparse tensors.
%
%   MRDIVIDE(A,B) is called for the syntax 'A / B' when A is a sparse
%   tensor and B is a scalar.
%
%   Example
%   X = sptenrand([4 3 2],5);
%   X / 3
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isscalar(B)
    newsubs = A.subs;
    newvals = A.vals / B;
    if B == 0
        nansubs = setdiff(allsubs(A),newsubs,'rows');
        newsubs = [newsubs; nansubs];
        newvals = [newvals; repmat(NaN,size(nansubs,1),1)];
    end
    C = sptensor(newsubs,newvals,A.size);
    return;
end

error('MRDIVIDE only supports the scalar case for sparse tensors');
