function z = isequal(x,y)
%ISEQUAL Compare spares tensors for equality.
%
%   ISEQUAL(A,B) compares the sparse tensors A and B for equality.
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%% Observations for sparse matrix case.
% The result of isequal(a,full(a)) is true!

%%
if ~isequal(x.size,y.size) 
    z = false;
elseif isa(x,'sptensor') && isa(y,'sptensor') 
    z = (nnz(x-y) == 0);
elseif isa(y,'tensor')
    z = isequal(full(x),y);
else
    z = false;
end
