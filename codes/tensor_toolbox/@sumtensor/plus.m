function t = plus( t, x )
%PLUS Plus for sumtensor.
%
%   T = T + A adds A to the sumtensor parts, where T is a sumtensor and 
%   A is any valid sumtensor component.
%
%   T = T + {A,B,C} adds A, B, and C to the sumtensor parts, assuming they
%   are valid sumtensor components.
%
%   Note that new parts are appended to the end of T, even if A + T is called
%
%   See also SUMTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

% The full license terms can be found in the file LICENSE.t


    
if iscell(x)
    t = sumtensor(t.part{:}, x{:});
else    
    t = sumtensor(t.part{:}, x);
end
