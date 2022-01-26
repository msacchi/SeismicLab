function b = permute(a,order)
%PERMUTE Permute dimensions for a ttensor.
%
%   Y = PERMUTE(X,ORDER) rearranges the dimensions of X so that they
%   are in the order specified by the vector ORDER. The tensor
%   produced has the same values of X but the order of the subscripts
%   needed to access any particular element are rearranged as
%   specified by ORDER.  
%
%   See also TTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


N = ndims(a);

if ~isequal(1:N,sort(order))
  error('Invalid permuation');
end

newcore = permute(a.core,order);
newu = a.u(order);
b = ttensor(newcore,newu);




