function Y = squeeze(X)
%SQUEEZE Remove singleton dimensions from a sparse tensor.
%
%   Y = SQUEEZE(X) returns a sparse tensor Y with the same elements as
%   X but with all the singleton dimensions removed.  A singleton
%   is a dimension such that size(X,dim)==1.  
%
%   If X has *only* singleton dimensions, then Y is a scalar.
%
%   Examples
%   squeeze( sptenrand([2,1,3],0.5) ) %<-- returns a 2-by-3 sptensor
%   squeeze( sptensor([1 1 1],1,[1 1 1]) ) %<-- returns a scalar
%   See also SPTENSOR.
% 
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if all(X.size > 1)
  % No singleton dimensions to squeeze
  Y = X;
else
  idx = find(X.size > 1);
  if numel(idx) == 0
    % Scalar case - only singleton dimensions
    Y = X.vals;
  else
    siz = X.size(idx);
    if isempty(X.vals)
        Y = sptensor([],[],siz);
    else
        Y = sptensor(X.subs(:,idx), X.vals, siz);
    end
  end
end

return;
