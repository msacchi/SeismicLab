function m = size(t,idx)
%SIZE Size of ktensor.
%  
%   D = SIZE(T) returns the size of the tensor. 
%
%   I = SIZE(T,DIM) returns the size of the dimension specified by
%   the scalar DIM.
%
%   See also KTENSOR, KTENSOR/NDIMS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isempty(t.lambda)
    m = [];
end

if exist('idx','var')
    m = size(t.u{idx}, 1);
else
    for i = 1 : ndims(t)
	m(i) = size(t.u{i}, 1);
    end
end
