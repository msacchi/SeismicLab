function t = full(t)
%FULL Convert a symktensor to a symtensor.
%
%   T = FULL(C) converts a symktensor to a symtensor.
%
%   Examples
%   X = symktensor([3; 2], ones(4,2));
%   Y = full(A) %<-- equivalent dense tensor
%
%   See also SYMKTENSOR, TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


n = size(t,1);
m = ndims(t);
sz = nchoosek(m+n-1,m);

I = zeros(sz,m);
for loc = 1:sz
    if loc == 1
        I(loc,:) = ones(1,m);
    else
        I(loc,:) = I(loc-1,:);
        j = m;
        while (I(loc,j) == n)
            j = j - 1;
        end
        I(loc,j:m) = I(loc,j)+1;
    end
end

tnew = symtensor(@ones,m,n);
vals = entry(t,I);
tnew(I) = vals;
t = tnew;
