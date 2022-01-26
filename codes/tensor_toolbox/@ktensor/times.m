function C = times(A,B)
%TIMES Element-wise multiplication for ktensor.
%
%   TIMES(A,B) denotes element-by-element multiplication (only supports the
%   second argument being a tensor or sptensor). 
% 
%   C = TIMES(A,B) is called for the syntax 'A .* B'. Either A or B must be
%   a tensor or sptensor.
%
%   See also KTENSOR, SPTENSOR/TIMES, TENSOR/TIMES.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if ~isequal(size(A),size(B))
    error('Must be two tensors of the same size');
end

switch class(B)
    case {'sptensor','tensor'}
        % Call back to sptensor version.
        C = times(B,A);
        return;
    otherwise
        error('Invalid second argument for ktensor/times');
end
