function Y = ttm(X,V,varargin)
%TTM Tensor times matrix.
%
%   Y = TTM(X,A,N) computes the n-mode product of tensor X with a
%   matrix A; i.e., X x_N A.  The integer N specifies the dimension
%   (or mode) of X along which A should be multiplied.  If size(A) =
%   [J,I], then X must have size(X,N) = I.  The result will be the
%   same order and size as X except that size(Y,N) = J.
%
%   Y = TTM(X,{A,B,C,...}) computes the n-mode product of tensor X
%   with a sequence of matrices in the cell array.  The n-mode
%   products are computed sequentially along all dimensions (or modes)
%   of X. The cell array contains ndims(X) matrices.
%
%   Y = TTM(X,{A,B,C,...},DIMS) computes the sequence tensor-matrix
%   products along the dimensions specified by DIMS.
%
%   Y = TTM(...,'t') performs the same computations as above except
%   the matrices are transposed.
%
%   Examples
%   X = tensor(rand(5,3,4,2));
%   A = rand(4,5); B = rand(4,3); C = rand(3,4); D = rand(3,2);
%   Y = ttm(X, A, 1)         %<-- computes X times A in mode-1
%   Y = ttm(X, {A,B,C,D}, 1) %<-- same as above
%   Y = ttm(X, A', 1, 't')   %<-- same as above
%   Y = ttm(X, {A,B,C,D}, [1 2 3 4]) %<-- 4-way multiply
%   Y = ttm(X, {D,C,B,A}, [4 3 2 1]) %<-- same as above
%   Y = ttm(X, {A,B,C,D})            %<-- same as above
%   Y = ttm(X, {A',B',C',D'}, 't')   %<-- same as above
%   Y = ttm(X, {C,D}, [3 4])     %<-- X times C in mode-3 & D in mode-4
%   Y = ttm(X, {A,B,C,D}, [3 4]) %<-- same as above
%   Y = ttm(X, {A,B,D}, [1 2 4])   %<-- 3-way multiply
%   Y = ttm(X, {A,B,C,D}, [1 2 4]) %<-- same as above
%   Y = ttm(X, {A,B,D}, -3)        %<-- same as above
%   Y = ttm(X, {A,B,C,D}, -3)      %<-- same as above
%
%   See also TENSOR, TENSOR/TTT, TENSOR/TTV.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>




%% Check the number of arguments
if (nargin < 2)
    error('TTM requires at least two arguments.');
end

%% Create 'n' and 'tflag' arguments from varargin
n = 1:ndims(X);
tflag = '';
ver = 0;
if numel(varargin) == 1
    if ischar(varargin{1})
        tflag = varargin{1};
    else
        n = varargin{1};
    end
elseif numel(varargin) == 2
    n = varargin{1};
    tflag = varargin{2};
elseif numel(varargin) == 3
    n = varargin{1};
    tflag = varargin{2};
    ver = varargin{3};
end

%% Handle cell array
if iscell(V)   
    % Copy n into dims
    dims = n;
    % Check that the dimensions are valid
    [dims,vidx] = tt_dimscheck(dims,ndims(X),numel(V));
    % Calculate individual products
    Y = ttm(X, V{vidx(1)}, dims(1), tflag);
    for k = 2 : numel(dims)
        Y = ttm(Y, V{vidx(k)}, dims(k), tflag);
    end
    % All done
    return;
end

%% Check the second argument
if ~ismatrix(V)
    error('tensor/ttm: 2nd argument must be a matrix.');
end

%% Check n
if (numel(n) ~= 1 || (n < 0) || n > ndims(X))
    error('Dimension N must be between 1 and NDIMS(X).');
end

%% COMPUTE SINGLE N-MODE PRODUCT 

N = ndims(X);
sz = size(X);

if ver == 0  %old verion
    order = [n,1:n-1,n+1:N];
    newdata = double(permute(X,order));
    newdata = reshape(newdata,sz(n),prod(sz([1:n-1,n+1:N])));
    if tflag == 't'
        newdata = V' * newdata;
        p = size(V,2);
    else
        newdata = V*newdata;
        p = size(V,1);
    end
    newsz = [p,sz(1:n-1),sz(n+1:N)];
    Y = tensor(newdata,newsz);
    Y = ipermute(Y,order);
else % new version
    
    if tflag == 't'
        p = size(V, 2);
    else
       p = size(V, 1);
    end

    
    if n == 1
        A = reshape(X.data, sz(n), []);
        if tflag == 't'
            B = V' * A;
        else
            B = V * A;
        end
    elseif n == N
        At = reshape(X.data, [], sz(n));
        if tflag == 't'
            B = At * V;
        else
            B = At * V';
        end
    else
        nblocks = prod(sz(n+1:N));
        ncols = prod(sz(1:n-1)); % Per block
        nAk = sz(n) * ncols; % Number entries in each block of A
        nBk = p * ncols;
        B = zeros(p * nblocks * ncols, 1);
        for k = 1 : nblocks
            % Extract k-th sub-block of A (in row-major orientation)
            Akt = reshape(X.data((k-1) * nAk + 1: k * nAk), ncols, sz(n));
            if tflag == 't'
                Bkt = Akt * V;
            else
                Bkt = Akt * V';
            end
            B((k-1) * nBk + 1: k * nBk) = Bkt(:);        
        end
        
    end
    newsz = sz;
    newsz(n) = p;
    Y = tensor(B, newsz);
   
   
end

return;

end
