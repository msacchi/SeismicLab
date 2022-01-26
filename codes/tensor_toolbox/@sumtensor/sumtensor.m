%SUMTENSOR Class for implicit sum of other tensors.
%
%SUMTENSOR Methods:
%   disp      - Command window display of a sumtensor.
%   display   - Command window display of a sumtensor.
%   double    - Convert sumtensor to double array.
%   full      - Convert a sumtensor to a (dense) tensor.
%   innerprod - Efficient inner product with a sumtensor.
%   isscalar  - False for sumtensors.
%   mttkrp    - Matricized tensor times Khatri-Rao product for sumtensor.
%   ndims     - Return the number of dimensions for a sumtensor.
%   norm      - Frobenius norm of a sumtensor.
%   plus      - Plus for sumtensor.
%   size      - Size of a sumtensor.
%   subsref   - Subscript reference for sumtensor.
%   sumtensor - Tensor stored as sum of tensors.
%   ttv       - Tensor times vector for sumtensor.
%   uminus    - Unary minus for sumtensor.
%   uplus     - Unary plus for sumtensor.
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','sumtensor_doc.html')))">Documentation page for sum of tensors class</a>
%
%   See also TENSOR_TOOLBOX
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

function t = sumtensor(varargin)
%SUMTENSOR Tensor stored as sum of tensors.
%
%   The SUMTENSOR class is limited to certain operations that easily
%   decompose as sums: INNERPROD, MTTKRP, TTV. Note that the NORM function
%   is not easily computed for a SUMTENSOR.
%
%   T = SUMTENSOR(T1,T2,...) creates a tensor that is the sum of its
%   constituent parts. The tensor is stored implicitly, i.e., each
%   component is retained. This may lead to storage and computation
%   efficiency. All input tensors must be the same size, but they can be
%   any type of tensor. 
%
%   T = SUMTENSOR(S) creates a SUMTENSOR by copying an existing
%   SUMTENSOR.
%
%   T = SUMTENSOR is the empty constructor.
%
%   Examples
%   T1 = tensor(rand(4,3,3));
%   T2 = sptensor([1 1 1; 3 1 2; 4 3 3], 1, [4,3,3]);
%   T = sumtensor(T1,T2); %<--A sumtensor with parts T1 and T2
%
%   See also TENSOR, SPTENSOR, TTENSOR, KTENSOR
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


% Empty constructor
if (nargin == 0)
    t.part = cell(0);
    t = class(t, 'sumtensor');
    return;
end

% Copy constructor
if (nargin == 1) && isa(varargin{1}, 'sumtensor')
    t.part = varargin{1}.part;
    t = class(t, 'sumtensor');
    return;
end

% Multiple arguments constructor
t.part = cell(nargin,1);
for i = 1:nargin
    cl = class(varargin{i});
    if ismember(cl,'double') % Convert an MDA
        varargin{i} = tensor(varargin{i});
    elseif ~ismember(cl, {'tensor','sptensor','ktensor','ttensor'})
        error('Inputs must be tensors. Symtensors are not supported.');
    end
    
    if (i > 1)
        if ~isequal(size(varargin{i}), size(varargin{1}))
            error('All inputs must be the same size.');
        end
    end
    
    t.part{i} = varargin{i};  
    
end
t = class(t, 'sumtensor');
return;
