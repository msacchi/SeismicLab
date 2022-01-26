function export_data(A, fname, varargin)
%EXPORT_DATA Export tensor-related data to a file.
%   
%   EXPORT_DATA(A,FNAME) exports object A to the file named FNAME in plain ASCII
%   text. Export currently supports exporting the following data types:
%
%      - tensor
%      - sptensor
%      - ktensor
%      - matrix
%
%   In the case of a tensor, the first three lines give details about the
%   tensor. The format for a 4 x 3 x 2 tensor is as follows...
%
%      tensor
%      3
%      4 3 2
%      <value of A(1,1,1)>
%      <value of A(2,1,1)>
%      <value of A(3,1,1)>
%      <value of A(4,1,1)>
%      <value of A(1,2,1)>
%      <value of A(2,2,1)>
%      <value of A(3,2,1)>
%      <value of A(4,2,1)>
%      ...
%      <value of A(2,3,2)>
%      <value of A(3,3,2)>
%      <value of A(4,3,2)>
%
%   In the case of an sptensor, the first four lines give details about the
%   sptensor. The format for a 4 x 3 x 2 sptensor with 10 nonzeros is as 
%   follows...
%
%      sptensor
%      3
%      4 3 2
%      10
%      i1 j1 k1 
%      i2 j2 k2 <value of A(i2,j2,k2)>
%      ...
%      i10 j10 k10 <value of A(i10,j10,k10)>
%
%   A matrix is formatted the same as a 2-way tensor except that the first
%   line says "matrix" rather than "tensor".
%
%   In the case of an ktensor, the first four lines give details about 
%   the ktensor. The format for a 3 x 4 x 5 ktensor with rank=2 is as 
%   follows, which contains the lambda values and factor matrices, U, 
%   (one per tensor dimension), ...
%
%      ktensor
%      3 
%      3 4 5 
%      2 
%      <value of lambda1> <value of lambda2>
%      matrix
%      2 
%      3 2 
%      <value of U{1}(1,1)> <value of U{1}(1,2) 
%      ...
%      <value of U{1}(3,1)> <value of U{1}(3,2) 
%      matrix
%      2 
%      4 2 
%      <value of U{2}(1,1)> <value of U{2}(1,2) 
%      ...
%      <value of U{2}(4,1)> <value of U{2}(4,2) 
%      matrix
%      2 
%      5 2 
%      <value of U{3}(1,1)> <value of U{3}(1,2) 
%      ...
%      <value of U{3}(5,1)> <value of U{3}(5,2) 
%
%   EXPORT_DATA(A,FNAME,'FMT_*',FORMAT) supports formatting of different
%   parts of the output using optional input arguments. FORMAT is a
%   character vector that follows the convention used in FPRINTF. The
%   different formatting strings are as follows:
%
%      FMT_DATA:   format string for data elements of all data [default: '%.16e']
%      FMT_LAMBDA: format string for lambda values of KTENSOR data [default: '%.16e']
%
%   Examples
%   X = sptensor([1 1 1;2 2 2; 4 3 2],[1;2;3],[4 3 2])
%   export_data(X,'X_default.sptensor') 
%   export_data(X,'X_int_data.sptensor','fmt_data','%d')
%   X1 = import_data('X_default.sptensor')
%   X2 = import_data('X_int_data.sptensor')
%   isequal(X1,X), isequal(X2,X), isequal(X1,X2)
%   K = ktensor([1; 2], 4*ones(4,2), 5*ones(5,2), 3*ones(3,2))
%   export_data(K,'K_default.ktensor') 
%   export_data(K,'K_int_data.ktensor','fmt_data','%d') 
%   export_data(K,'K_int_data_int_lambda.ktensor','fmt_data','%d','fmt_lambda','%d') 
%   K1 = import_data('K_default.ktensor')
%   K2 = import_data('K_int_data.ktensor')
%   K3 = import_data('K_int_data_int_lambda.ktensor')
%   isequal(K1,K), isequal(K2,K), isequal(K3,K)
%
%   See also TENSOR, SPTENSOR, KTENSOR, IMPORT_DATA, FPRINTF
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

%% Parse the inputs
p = inputParser;
addRequired(p,'A',@isobject);
addRequired(p,'fname',@ischar);
addParameter(p,'fmt_data', '%.16e', @ischar);
addParameter(p,'fmt_lambda', '%.16e', @ischar);
parse(p,A,fname,varargin{:});
A = p.Results.A;
fname = p.Results.fname;
fmt_data = p.Results.fmt_data;
fmt_lambda = p.Results.fmt_lambda;


%% Open file
fid = fopen(fname,'W');
if (fid == -1)
    error('Cannot open file %s',fname);
end


%% Export the object

if isa(A,'tensor')
    
    fprintf(fid, 'tensor\n');
    export_size(fid, size(A));
    export_array(fid, A.data, fmt_data);   
 
elseif isa(A,'sptensor')
    
    fprintf(fid, 'sptensor\n');
    export_sparse_size(fid, A);
    export_sparse_array(fid, A, fmt_data);   

elseif isa(A,'ktensor')
    
    fprintf(fid, 'ktensor\n');
    export_size(fid, size(A));
    export_rank(fid, A);
    export_lambda(fid, A.lambda, fmt_lambda);   
    for n = 1:length(size(A))
        fprintf(fid, 'matrix\n');
        export_size(fid, size(A.U{n}));
        export_factor(fid, A.U{n}, fmt_data);
    end

elseif isnumeric(A) && ndims(A) == 2        

    fprintf(fid, 'matrix\n');
    export_size(fid, size(A));
    export_array(fid, A, fmt_data);    
    
else   
    
    error('Invalid data type for export');    
    
end


%% Close file
fclose(fid);

function export_size(fid, sz)
% Export the size of something to a file
fprintf(fid, '%d \n', length(sz)); % # of dimensions on one line
fprintf(fid, '%d ', sz); % # size of each dimensions on the next line
fprintf(fid, '\n');

function export_rank(fid, data)
% Export the rank of a ktensor to a file
fprintf(fid, '%d \n', length(data.lambda)); % ktensor rank on one line

function export_lambda(fid, data, fmt_data)
% Export dense data that supports numel and linear indexing
fprintf(fid, append(fmt_data,' '), data);
fprintf(fid, '\n');

function export_array(fid, data, fmt_data)
% Export dense data that supports numel and linear indexing
for i = 1:numel(data)
    fprintf(fid, append(fmt_data,'\n'), data(i));
end

function export_factor(fid, data, fmt_data)
% Export dense data that supports numel and linear indexing
for i = 1:size(data,1)
    fprintf(fid, append(fmt_data,' '), data(i,:));
    fprintf(fid, '\n');
end

function export_sparse_size(fid, A)
% Export the size of something to a file
fprintf(fid, '%d \n', length(size(A))); % # of dimensions on one line
fprintf(fid, '%d ', size(A)); % # size of each dimensions on the next line
fprintf(fid, '\n');
fprintf(fid, '%d \n', nnz(A)); % # number of nonzeros on the next line

function export_sparse_array(fid, A, fmt_data)
% Export sparse array data in coordinate format
fmt_str = append(repmat('%d ',1,length(size(A))),fmt_data);
data = [A.subs A.vals];
data_str = compose(fmt_str,data);
fprintf(fid,'%s\n',data_str{:}); 
