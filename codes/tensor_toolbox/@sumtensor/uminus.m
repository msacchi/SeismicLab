function t = uminus(t)
%UMINUS Unary minus for sumtensor.
%
%   See also SUMTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


for i = 1:lenght(t.part)
    t.part{i} = -t.part{i};
end


