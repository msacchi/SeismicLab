function res = innerprod(X,Y)
%INNERPROD Efficient inner product with a ttensor.
%
%   R = INNERPROD(X,Y) efficiently computes the inner product between
%   two tensors X and Y. This inner product is the standard inner product,
%   if the tensors were treated as vectors. How to do this most efficiently 
%   depends on the tensor Y.
%
%   See also TENSOR/INNERPROD, TTENSOR, KTENSOR/INNERPROD
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


% X is a ttensor
switch class(Y)
    
    case {'ttensor'}
        if ~isequal(size(X),size(Y))
            error('X and Y must be the same size.');
        end
        if prod(size(X.core)) > prod(size(Y.core))
            % Reverse argument and call this function again so that the
            % tensor with the smaller core is the first argument. 
            res = innerprod(Y,X);
            return
        end
        W = cell(ndims(X),1);
        for n = 1:ndims(X)
            W{n} = X.u{n}'*Y.u{n};
        end
        J = ttm(Y.core, W);
        res = innerprod(X.core,J);
        return

    case {'tensor','sptensor'}
        if ~isequal(size(X),size(Y))
            error('X and Y must be the same size.');
        end
        if (prod(size(X)) < prod(size(X.core)))
            Z = full(X);
            res = innerprod(Z,Y);
            return;
        end
        Z = ttm(Y,X.u,'t');
        res = innerprod(Z, X.core);
        return

    case {'ktensor'}
        % Reverse arguments to call ktensor implementation
        res = innerprod(Y,X);
        return
        
    otherwise
        disp(['Inner product not available for class ' class(Y)]);
        return
end


