function sz = size(t,idx)
%SIZE Dimensions of a symmetric tensor.
%  
%   D = SIZE(T) returns the size of a symtensor T. 
%
%   I = SIZE(T,DIM) returns the size of symtensor T in the dimension 
%   specified by the scalar DIM.
%
%   See also SYMTENSOR, SYMTENSOR/NDIMS, SIZE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if (t.m == 0)
    sz = [];
elseif exist('idx','var')
    if 1<=idx && idx<=t.m  %Bounds check
        sz = t.n;
    else
        error('Index exceeds tensor dimensions');
    end
else
    sz = t.n * ones(1,t.m);  
end
