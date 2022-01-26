function G = tenrandblk(bsz, bns, verbose)
%TENRANDBLK Generate nearly block diagonal tensor.
%
%   G = TENRANDBLK(BSZ,BNS) creates a tensor G that is block 'diagonal'
%   plus noise. The first argument specifies the size of each block as a
%   row. The number of rows of BSZ is the number of blocks and the order of
%   the tensor is the number of columns. The blocks need not be square. The
%   size of G is equal to sum(BSZ,1). The second argument specifies the
%   squared norm of each block. The values must be strictly decreasing and
%   sum to less than one. The squared norm of the offdiagonal parts of G is
%   1-sum(BNS).  
%
%   See <a
%   href="matlab:p=what('tensor_toolbox');web(strcat('file:',p.path,'/doc/html/T5_hosvd_algorithm_doc.html'))">example usage</a>.
%
%   See also CREATE_PROBLEM
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

%% Check inputs

% --- Check sizes ---

% Must be a matrix
if ~ismatrix(bsz)
    error('Block sizes must be a matrix');
end

if any(bsz <= 0)
    error('Block sizes must be positive');
end

if any(bsz ~= round(bsz))
    error('Block sizes must be integers');
end

% Find the end of each block
bend = cumsum(bsz,1);

% Extract size: D = # dimenions, L = # levels
[D,L] = size(bsz);

% Final size
gsz = bend(end,:);

% --- Check errors ---

% Check errors are okay
if sum(bns) > 1
    error('sum of all block squared error norms must be <= 1');
end

for i = 2:L
    if bns(i) > bns(i-1)
        error('Squared norms must be strictly decreasing');
    end
end

if ~exist('verbose','var')
    verbose = false;
end

%% Create tensor

% Figure out norm of off-block-diagonal
dltnrmsqr = 1 - sum(bns);

% Create pattern for off-block-diagonal to be modified as we go
dltpattern = ones(gsz);

% Create tensor to fill in
G = tensor(@zeros,gsz);

% Create random entries to use 
Grnd = tensor(@(sz) sign(randn(sz)) .* (0.1*rand(sz)+0.9), gsz);

% Loop through and create blocks
for i = 1:L
    
    % Figure out ith block pattern
    blkrange = cell(D,1);
    for k = 1:D
        if i == 1
            blkrange{k} = 1:bend(i,k);
        else
            blkrange{k} = bend(i-1,k)+1:bend(i,k);
        end
    end
    
    % Create pattern that has ones for the block
    pattern = zeros(gsz);
    pattern(blkrange{:}) = 1;
    
    % Zero out block in the off-diagonal pattern
    dltpattern(blkrange{:}) = 0;
      
    % Randomly fill delta-pattern and rescale
    block = Grnd .* pattern;
    sse = collapse(block.^2);
    block = sqrt(bns(i)/sse) .* block;

    % Add to main tensor
    G = G + block;

    % Verbose output
    if verbose
        fprintf('Created block of size %s with norm (%f)^2=%f\n', tt_size2str(bsz(i,:)), norm(block), norm(block)^2);
    end
end

if dltnrmsqr > 0
    % Final pattern
    block = Grnd .* dltpattern;
    sse = collapse(block.^2);
    block = sqrt(dltnrmsqr/sse) .* block;
    G = G + block;
    if verbose
        fprintf('Created tensor of size %s with off-block-diaognal norm (%f)^2=%f\n', tt_size2str(gsz), norm(block), norm(block)^2);
    end
end
