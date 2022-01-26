%SPTENSOR Class for sparse tensors.
%
%SPTENSOR Methods:
%   and       - Logical AND (&) for sptensors.
%   collapse  - Collapse sparse tensor along specified dimensions.
%   contract  - Contract sparse tensor along two dimensions (array trace).
%   disp      - Command window display of a sparse tensor.
%   display   - Command window display of a sparse tensor.
%   divide    - Divide an SPTENSOR by a nonnegative KTENSOR.
%   double    - Converts a sparse tensor to a dense multidimensional array.
%   elemfun   - Manipulate the nonzero elements of a sparse tensor.
%   end       - Last index of indexing expression for sparse tensor.
%   eq        - Equal (==) for sptensors.
%   find      - Find subscripts of nonzero elements in a sparse tensor.
%   full      - Convert a sparse tensor to a (dense) tensor.
%   ge        - Greater than or equal for sptensors.
%   gt        - Greater than for sptensors.
%   innerprod - Efficient inner product with a sparse tensor.
%   isequal   - Compare spares tensors for equality.
%   isscalar  - False for sptensors.
%   ldivide   - Array right division for sparse tensors.
%   le        - Less than or equal for sptensors.
%   lt        - Less than for sptensors.
%   mask      - Extract values as specified by a mask tensor.
%   minus     - Binary subtraction for sparse tensors. 
%   mldivide  - Slash left division for sparse tensors.
%   mrdivide  - Slash right division for sparse tensors.
%   mtimes    - sptensor-scalar multiplication.
%   mttkrp    - Matricized tensor times Khatri-Rao product for sparse tensor.
%   ndims     - Number of dimensions of a sparse tensor.
%   ne        - Not equal (~=) for sptensors.
%   nnz       - Number of nonzeros in sparse tensor.
%   norm      - Frobenius norm of a sparse tensor.
%   not       - Logical NOT (~) for sptensors.
%   nvecs     - Compute the leading mode-n vectors for a sparse tensor.
%   ones      - Replace nonzero elements of sparse tensor with ones.
%   or        - Logical OR (|) for sptensors.
%   permute   - Rearrange the dimensions of a sparse tensor.
%   plus      - Binary addition for sparse tensors. 
%   rdivide   - Array right division for sparse tensors.
%   reshape   - Reshape sparse tensor.
%   scale     - Scale along specified dimensions for sparse tensors.
%   size      - Sparse tensor dimensions.
%   spmatrix  - Converts a two-way sparse tensor to sparse matrix.
%   spones    - Replace nonzero sparse tensor elements with ones.
%   sptensor  - Create a sparse tensor.
%   squeeze   - Remove singleton dimensions from a sparse tensor.
%   subsasgn  - Subscripted assignment for sparse tensor.
%   subsref   - Subscripted reference for a sparse tensor.
%   times     - Array multiplication for sparse tensors.
%   ttm       - Sparse tensor times matrix.
%   ttt       - Sparse tensor times sparse tensor.
%   ttv       - Sparse tensor times vector.
%   uminus    - Unary minus (-) for sptensor.
%   uplus     - Unary plus (+) for sptensor.
%   xor       - Logical XOR for sptensors.
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','sptensor_doc.html')))">Documentation page for Sparse Tensor Class</a>
%
%   See also TENSOR_TOOLBOX
%
%   How to cite the sptensor class:
%   * B.W. Bader and T.G. Kolda. Efficient MATLAB Computations with Sparse
%     and Factored Tensors, SIAM J. Scientific Computing, 30:205-231, 2007,
%     http://dx.doi.org/10.1137/060676489. 
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

function t = sptensor(varargin)
%SPTENSOR Create a sparse tensor.
%
%   X = SPTENSOR(SUBS, VALS, SZ, FUN) uses the rows of SUBS and VALS
%   to generate a sparse tensor X of size SZ = [m1 m2 ... mn]. SUBS is
%   an p x n array specifying the subscripts of the values to be
%   inserted into S. The k-th row of SUBS specifies the subscripts for
%   the k-th value in VALS. The values are accumulated at repeated
%   subscripts using the function FUN, which is specified by a
%   function handle.
%
%   There are several simplifications of this four argument call.
%
%   X = SPTENSOR(SUBS,VALS,SZ) uses FUN=@SUM.
%
%   X = SPTENSOR(SUBS,VALS) uses SM = max(SUBS,[],1).
%
%   X = SPTENSOR(SZ) abbreviates X = SPTENSOR([],[],SZ).
%
%   X = SPTENSOR(Y) copies/converts Y if it is an sptensor, an sptenmat, or
%   a dense tensor or MDA (the zeros are squeezed out), an sptensor3, or a
%   sparse matrix. Note that a row-vector, integer MDA is interpreted as a
%   size (see previous constructor).
%
%   X = SPTENSOR is the empty constructor.
%
%   X = SPTENSOR(FH,SZ,NZ) creates a random sparse tensor of the specified
%   size with NZ nonzeros (this can be an explit value or a proportion).
%   The function handle FH is used to create the nonzeros.
%
%   The argument VALS may be scalar, which is expanded to be the
%   same length as SUBS, i.e., it is equivalent to VALS*(p,1).
%
%   Examples
%   subs = [1 1 1; 1 1 3; 2 2 2; 4 4 4; 1 1 1; 1 1 1]
%   vals = [0.5; 1.5; 2.5; 3.5; 4.5; 5.5]
%   siz = [4 4 4];
%   X = sptensor(subs,vals,siz) %<-- sparse 4x4x4, repeats summed
%   X = sptensor(subs,1,siz) %<-- scalar 2nd argument
%   X = sptensor(subs,vals,siz,@max) %<-- max for accumulation
%   myfun = @(x) sum(x) / 3;
%   X = sptensor(subs,vals,siz,myfun) %<-- custom accumulation
%
%   See also SPTENSOR, SPTENRAND.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



% EMPTY Constructor
if (nargin == 0) || ((nargin == 1) && isempty(varargin{1}))
    t.subs = [];
    t.vals = [];
    t.size = [];
    t = class(t,'sptensor');
    return;
end

% SINGLE ARGUMENT
if (nargin == 1)

    source = varargin{1};

    switch(class(source))

        % COPY CONSTRUCTOR
        case 'sptensor',                
            t.subs = source.subs;
            t.vals = source.vals;
            t.size = source.size;
            t = class(t, 'sptensor');
            return;
            
        % CONVERT SPTENMAT
        case 'sptenmat',                

            % Extract the tensor size and order
            siz = source.tsize;
            
            if isempty(source.subs) %There are no nonzero terms                
                subs = [];            
            else % Convert the 2d-subscipts into nd-subscripts                                
                if ~isempty(source.rdims)
                    subs(:,source.rdims) = ... 
                        tt_ind2sub(siz(source.rdims),source.subs(:,1));
                end
                if ~isempty(source.cdims)
                    subs(:,source.cdims) = ...
                        tt_ind2sub(siz(source.cdims),source.subs(:,2));
                end
            end
            % Copy the values (which do not need to be modified)
            vals = source.vals;

            % Store everything
            t.subs = subs;
            t.vals = vals;
            t.size = siz;
            t = class(t, 'sptensor');
            return;

        % CONVERT TENSOR
        case 'tensor',                  
            [subs,vals] = find(source); 
            t.subs = subs;
            t.vals = vals;    
            t.size = size(source);
            t = class(t, 'sptensor');
            return;

        % CONVERT SPTENSOR3
        case 'sptensor3',
            K = size(source,3);
            [I,J] = size(source{1});
            nz = nnz(K);  
            bigsubs = [];
            bigvals = [];
            for k = 1:K
                [subs,vals] = find(source{k});
                if isempty(bigsubs)                 
                    bigsubs = [subs, k*ones(size(subs,1),1)];
                    bigvals = [vals];
                else
                    bigsubs = [bigsubs; subs, k*ones(size(subs,1),1)];
                    bigvals = [bigvals; vals];
                end
            end
            t.subs = bigsubs;
            t.vals = bigvals;
            t.size = [ I J K ];
            t = class(t,'sptensor');               
            return;
            
        % SPARSE MATRIX, SIZE, or MDA
        case {'numeric','logical','double'},

            % Case 1: SPARSE MATRIX
            if issparse(source)
                [i,j,s] = find(source);
                siz = size(source);
                t.subs = [i,j];
                t.vals = s;
                t.size = siz;
                t = class(t,'sptensor');
                return;
            end

            % Case 2: SPECIFYING THE SIZE
            if tt_sizecheck(source)
                t.subs = [];
                t.vals = [];
                t.size = source;
                t = class(t, 'sptensor');
                return;
            end

            % Case 3: An MDA
            t = sptensor(tensor(source));
            return;

    end % switch

end % nargin == 1

% SPECIAL CASE for INTERACTION WITH MEX FILES OR DIRECT CREATION OF
% SPTENSOR WITHOUT ANY SORTING OR OTHER STANDARD CHECKS
if (nargin == 4) && (isnumeric(varargin{4})) && (varargin{4} == 0)

    % Store everything
    t.subs = varargin{1};
    t.vals = varargin{2};
    t.size = varargin{3};

    % Create the tensor
    t = class(t, 'sptensor');

    return;

end

% RANDOM TENSOR
if (nargin == 3) && isa(varargin{1},'function_handle')
    fh = varargin{1};
    sz = varargin{2};
    nz = varargin{3};
    
    if (nz < 0) || (nz >= prod(sz))
        error('Requested number of nonzeros must be positive and less than the total size')
    elseif (nz < 1)
        nz = ceil(prod(sz) * nz);
    else
        nz = floor(nz);
    end
    
    % Keep iterating until we find enough unique nonzeros or we give up
    subs = [];
    cnt = 0;
    while (size(subs,1) < nz) && (cnt < 10)
        subs = ceil( rand(nz, size(sz,2)) * diag(sz) );
        subs = unique(subs, 'rows');
        cnt = cnt + 1;
    end
    
    nz = min(nz, size(subs,1));
    subs = subs(1:nz,:);
    vals = fh(nz,1);
    
    % Store everything
    t.subs = subs;
    t.vals = vals;
    t.size = sz;

    % Create the tensor
    t = class(t, 'sptensor');
    return;
end

% CONVERT A SET OF INPUTS
if (nargin == 2) || (nargin == 3) || (nargin == 4)

    % Extract the subscripts and values
    subs = varargin{1};   
    vals = varargin{2};
    
    tt_subscheck(subs);
    tt_valscheck(vals);
    if ~isempty(vals) && (numel(vals) ~= 1) && (size(vals,1) ~= size(subs,1))
        error('Number of subscripts and values must be equal');
    end

    % Extract the size
    if (nargin > 2)
        siz = varargin{3};
        tt_sizecheck(siz);
    else
        siz = max(subs,[],1);
    end

    % Check for wrong input
    if size(subs,2) > size(siz,2)
        error('More subscripts than specified by size')
    end

    % Check for subscripts out of range
    for j = 1:numel(siz)
        if ~isempty(subs) && max(subs(:,j)) > siz(j)
            error('Subscript exceeds sptensor size')
        end
    end

    % Extract the 'combiner' function handle
    if (nargin == 4)
        fun = varargin{4};
    else
        fun = @sum;
    end
    
    if isempty(subs)
        newsubs = [];
        newvals = [];
    else
        % Identify only the unique indices
        [newsubs,junk,loc] = unique(subs,'rows');

        % Sum the corresponding values
        newvals = accumarray(loc,vals,[size(newsubs,1) 1],fun);
    end

    % Find the nonzero indices of the new values
    nzidx = find(newvals);
    newsubs = newsubs(nzidx,:);
    newvals = newvals(nzidx);

    % Store everything
    t.subs = newsubs;
    t.vals = newvals;
    t.size = siz;

    % Create the tensor
    t = class(t, 'sptensor');

    return;
end

error('Unsupported use of function SPTENSOR.');

