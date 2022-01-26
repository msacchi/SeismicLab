function y = not(x)
%NOT Logical NOT (~) for sptensors.
%
%   ~X performs a logical not on the input tensor X. The result always
%   returned as a sparse tensor.
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%% Observations for sparse matrix case.
% The result of ~a is sparse.

%% Then compute those indicies that are not in x
subs = setdiff(allsubs(x),x.subs,'rows');

%% Assemble final result
y = sptensor(subs,true,x.size);
