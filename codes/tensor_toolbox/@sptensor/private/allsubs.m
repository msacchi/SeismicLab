function s = allsubs(x)
%ALLSUBS Generate all possible subscripts for a sparse tensor X.
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%% Generate all possible indicies

% Preallocate (discover any memory issues here!)
s = zeros(prod(x.size),ndims(x));

% Generate appropriately sized ones vectors.
o = cell(ndims(x),1);
for n = 1:ndims(x)
    o{n} = ones(size(x,n),1);
end

% Generate each column of the subscripts in turn
for n = 1:ndims(x)
    i = o;
    i{n} = (1:size(x,n))';
    s(:,n) = khatrirao(i); 
end
