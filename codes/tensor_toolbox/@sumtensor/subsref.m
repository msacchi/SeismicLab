function B = subsref(A, S) 
%SUBSREF Subscript reference for sumtensor.
%
%   T.part{i} returns the ith part of the sumtensor T
%
%   Examples
%   T1 = tensor(rand(3,3,3));
%   T2 = sptensor([1 1 1; 3 1 2; 1 1 3], 1, [3,3,3]);
%   T = sumtensor(T1,T2); 
%   T.part{2} %<--Returns the symmetric tensor
%
%   See also SUMTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


B = builtin('subsref', A, S);

