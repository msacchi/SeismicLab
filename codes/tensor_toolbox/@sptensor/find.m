function [subs,vals] = find(t)
%FIND Find subscripts of nonzero elements in a sparse tensor.
%
%   [SUBS,VALS] = FIND(T) returns the subscripts and corresponding
%   values of the nonzero elements of T.
%
%   Note that unlike the standard MATLAB find function for an array,
%   find does not return linear indices. Instead, it returns an M x N
%   array where M is the number of nonzero values and N = ndims(T).
%   Thus, I(k,:) specifies the subscript of value V(k).
%
%   See also SPTENSOR, FIND.
% 
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



subs = t.subs;
vals = t.vals;
