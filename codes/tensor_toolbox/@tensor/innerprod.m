function res = innerprod(X,Y)
%INNERPROD Efficient inner product with a tensor.
%
%   R = INNERPROD(X,Y) efficiently computes the inner product between
%   two tensors X and Y. If Y is a tensor, then the inner product is
%   computed directly.  Otherwise, the inner product method for
%   that type of tensor is called.
%
%   Examples
%   X = tensor(rand(3,4,3));
%   Y = tensor(randi(10,3,4,3));
%   innerprod(X,Y)
%
%   See also TENSOR, SPTENSOR/INNERPROD, KTENSOR/INNERPROD, TTENSOR/INNERPROD
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



% X is a tensor
switch class(Y)
 
  case {'tensor'}
    % No need for same size check because it is implicit in the inner
    % product below. 
    if ~isequal(size(X), size(Y))
        error('TTB:UnequalSize', 'X and Y must be the same size');
    end
    x = reshape(X.data, 1, numel(X.data));
    y = reshape(Y.data, numel(Y.data), 1);
    res = x*y;
    
  case {'sptensor','ktensor','ttensor'}
    % Reverse arguments to call specialized code
    res = innerprod(Y,X);  
 
  otherwise
    disp(['Inner product not available for class ' class(Y)]);

end

return;
