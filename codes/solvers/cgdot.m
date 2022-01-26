function [out] = cgdot(in1,in2);
%cgdot: This is a flexible dot product function for cgls. I use this one 
%       when I work with linear operators with input/output that are not
%       vectors.
%
%  IN   in1,in2: vectors, matrices, cubes etc of the same size
%
%  OUT  out:  inner product (real scalar)
%
%
%    
% Mauricio D Sacchi, 2016 


 temp  =   in1.*conj(in2);
 out = sum(temp(:));

return;


