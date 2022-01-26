function z = isequal(x,y)
%ISEQUAL for symmetric tensors.
%
%   ISEQUAL(A,B) compares the symmetric tensors A and B for equality.
%
%   See also SYMTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%%
if ~isequal(size(x),size(y)) 
    z = false;
elseif isa(x,'symtensor') && isa(y,'symtensor') 
    z = isequal(x.val,y.val);
else
    z = false;
end
