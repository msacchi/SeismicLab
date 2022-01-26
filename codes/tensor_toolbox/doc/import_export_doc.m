%% Importing and Exporting Tensor Data
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="working.html">Working with Tensors</a> 
% &#62;&#62; <a href="import_export_doc.html">Importing and Exporting Tensor Data</a>
% </p>
% </html>
%
% You can import and export tensor data in ASCII format for the following
% data types:
%
% <html>
%   <ul>
%     <li><tt>tensor</tt></li>
%     <li><tt>matrix</tt></li>
%     <li><tt>sptensor</tt></li>
%     <li><tt>ktensor</tt></li>
%   </ul>
% </html>

%% Tensor/Matrix Data
% In the case of a tensor, the first three lines give details about the
% tensor. The format for a 4 x 3 x 2 tensor is as follows...
%
%    tensor
%    3
%    4 3 2
%    <value of A(1,1,1)>
%    <value of A(2,1,1)>
%    <value of A(3,1,1)>
%    <value of A(4,1,1)>
%    <value of A(1,2,1)>
%    <value of A(2,2,1)>
%    <value of A(3,2,1)>
%    <value of A(4,2,1)>
%    ...
%    <value of A(2,3,2)>
%    <value of A(3,3,2)>
%    <value of A(4,3,2)>
%
% A matrix is formatted the same as a 2-way tensor except that the first
% line says |matrix| rather than |tensor|.
%
Y1 = tenrand([4 3 2]);
export_data(Y1,'Y.tensor');
Y2 = import_data('Y.tensor');
assert(isequal(Y1,Y2)) %<-- test if equal when reading tensor from file
type 'Y.tensor'

%% Sptensor Data
% In the case of an sptensor, the first four lines give details about the
% sptensor. The format for a 4 x 3 x 2 sptensor with 10 nonzeros is as 
% follows...
%
%    sptensor
%    3
%    4 3 2
%    3
%    i1 j1 k1 
%    i2 j2 k2 <value of A(i2,j2,k2)>
%    i3 j3 k3 <value of A(i3,j3,k3)>
%
X1 = sptensor([1 1 1;2 2 2; 4 3 2],rand(3,1),[4 3 2]);
export_data(X1,'X.sptensor');
X2 = import_data('X.sptensor');
assert(isequal(X1,X2)) %<-- test if equal when reading tensor from file
type 'X.sptensor'

%% Ktensor Data
% In the case of an ktensor, the first four lines give details about 
% the ktensor. The format for a 4 x 3 x 2 ktensor with rank=2 is as 
% follows, which contains the lambda values and factor matrices, U, 
% (one per tensor dimension), ...
%
%    ktensor
%    3 
%    4 3 2 
%    2 
%    <value of lambda1> <value of lambda2>
%    matrix
%    2 
%    4 2 
%    <value of U{1}(1,1)> <value of U{1}(1,2) 
%    ...
%    <value of U{1}(4,1)> <value of U{1}(4,2) 
%    matrix
%    2 
%    3 2 
%    <value of U{2}(1,1)> <value of U{2}(1,2) 
%    ...
%    <value of U{2}(3,1)> <value of U{2}(3,2) 
%    matrix
%    2 
%    2 2 
%    <value of U{3}(1,1)> <value of U{3}(1,2) 
%    <value of U{3}(2,1)> <value of U{3}(2,2) 
%
K1 = ktensor(rand(2,1), rand(4,2), rand(3,2), rand(2,2));
export_data(K1,'K.ktensor');
K2 = import_data('K.ktensor');
assert(isequal(K1,K2)) %<-- test if equal when reading tensor from file
type 'K.ktensor'

%% Formatting Output - Data Values
% When exporting tensor data your can specify the format of the data values
% using character vector that follows the convention used in |fprintf|. For
% example, this can be useful when exporting tensors that contain only
% integers, resulting in a file that is smaller in size. The optional
% |fmt_data| parameter to |export_data| specifies the format string, with
% |'%.16|' being the default value.
%
X3 = sptensor([1 1 1;2 2 2; 4 3 2],[1;2;3],[4 3 2]);
export_data(X3,'X_fmt_data.sptensor','fmt_data','%d');
X4 = import_data('X_fmt_data.sptensor');
assert(isequal(X3,X4)) %<-- test if equal when reading tensor from file
type 'X_fmt_data.sptensor'

%% Formatting Output - Ktensor Lambda Values
% When exporting |ktensor| data your can specify the format of the |lambda|
% values using character vector that follows the convention used in
% |fprint|. For example, this can be useful when exporting tensors that
% contain only integers, resulting in a file that is smaller in size. The
% optional |fmt_lambda| parameter to |export_data| specifies the format
% string, with |'%.16'| being the default value.
%
K3 = ktensor([1; 2], 4*ones(4,2), 3*ones(3,2), 2*ones(2,2));
export_data(K3,'K_fmt_data_lambda.ktensor','fmt_data','%d','fmt_lambda','%d');
K4 = import_data('K_fmt_data_lambda.ktensor');
assert(isequal(K3,K4)) %<-- test if equal when reading tensor from file
type 'K_fmt_data_lambda.ktensor'
