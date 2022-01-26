function z = isequal(x,y)
%ISEQUAL for tensors.
%
%   ISEQUAL(A,B) compares the tensors A and B for equality.
%
%   See also TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%%
if ~isequal(x.size,y.size) 
    z = false;
elseif isa(x,'tensor') && isa(y,'tensor') 
    z = isequal(x.data,y.data);
elseif isa(y,'sptensor')
    z = isequal(x,full(y));
else
    z = false;
end
