function Z = minus(X,Y)
%MINUS Binary subtraction (-) for tenmat.
%
%   See also TENMAT, TENMAT/TENMATFUN.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



fun = @minus;

% One argument is a scalar
if ((prod(size(X)) == 1 || prod(size(Y)) == 1))
  if (prod(size(Y)) == 1) && isa(X,'tenmat')
    Z = X;
    Z.data = fun(Z.data,Y);
  else
    Z = Y;
    Z.data = fun(X,Z.data);
  end
  return;
end


% Both arguments are tenmats
Z = tenmat(Y);
if ~(isequal(size(Y),size(Z)))
  error('Tenmat size mismatch.')
end
Z.data = fun(X.data,Z.data);
