function b = permute(a,order)
%PERMUTE Permute dimensions of a ktensor.
%
%   B = PERMUTE(A,ORDER) rearranges the dimensions of A so that they
%   are in the order specified by the vector ORDER. The output is a ktensor
%   with components rearranged as specified by ORDER. The corresponding 
%   tensor has the same components as A but the order of the subscripts
%   needed to access any particular element is rearranged as specified by 
%   ORDER.
%
%   See also KTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



N = ndims(a);

if ~isequal(1:N,sort(order))
  error('Invalid permuation');
end

b = ktensor(a.lambda, a.u(order));




