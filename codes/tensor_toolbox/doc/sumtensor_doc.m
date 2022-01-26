%% Sum of Structured Tensors
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="tensor_types.html">Tensor Types</a> 
% &#62;&#62; <a href="sumtensor_doc.html">Sum of Structured Tensors</a>
% </p>
% </html>
%
% When certain operations are performed on a tensor which is formed as a 
% sum of tensors, it can be beneficial to avoid explicitly forming the sum.
% For example, if a tensor is formed as a sum of a low rank tensor and a 
% sparse tensor, the structure of the summands can make storage, decomposition and
% operations with other tensors significantly more efficient. The tensor
% toolbox supports a |sumtensor| object designed to exploit this structure.
% Here we explain the basics of defining and using sumtensors.
%% Creating sumtensors
% A sumtensor T can only be delared as a sum of same-sized tensors |T1,
% T2,...,TN|. The summand tensors are stored in a cell array, which define
% the "parts" of the sumtensor. The parts of a sumtensor can be (generic) 
% tensors (as |tensor|), sparse tensors (as |sptensor|), Kruskal tensors 
% (as |ktensor|), or Tucker tensors (as |ttensor|). An example of the use
% of the |sumtensor| constructor follows.
T1 = tensor(ones(3,3,3)); %<--A tensor
T2 = sptensor([1 1 1; 2 2 2; 3 3 2; 2 1 1], 1, [3,3,3]); %<--A sparse tensor

T = sumtensor(T1,T2)

%% An Large-Scale Example
% For large-scale problems, the |sumtensor| class may make the difference
% as to whether or not a tensor can be stored in memory. Consider the
% following example, where $\mathcal{T}$ is of size $1000 x 1000 x 1000$, 
% formed from the sum of a |ktensor| and an |sptensor|.
X1 = rand(500, 3); %Generating some factor matrices
X2 = rand(500, 3); 
X3 = rand(500, 3);
K = ktensor([1; 1; 1], X1, X2, X3);
S = sptenrand([500, 500, 500], 1e-100);

ST = sumtensor(K,S); %<-- Declare the sumtensor
TT = full(ST); %<-- Form the sum of the tensors explicitly

whos ST TT %<--Output the storage information for these variables

%%
% The difference in memory between the full and sumtensor is a factor of 10^5!
% Hence we prefer to use the sumtensor object whenever possible.
%% Further examples of the sumtensor constructer
% The sumtensor constructor can declare an empty sumtensor object, having 
% no parts, as follows
P = sumtensor()
%%
% |sumtensor| also supports use as a copy constructor.
S = sumtensor(P)
%% Use ndims and size for the dimensions of a sumtensor
% For a given sumtensor, |ndims| returns the number of modes of a sumtensor.
% Similarly, |size| returns a size array of the sumtensor.
ndims(T)
size(T)
%% Use full to convert a sumtensor to a "generic" tensor
% The |full| function can be used to convert a sumtensor to a generic tensor. Note that
% for large-scale tensors, this can a large amount of memory because each part of
% the sumtensor will be expanded and then summed.
full(T)
%% Use double to convert a sumtensor to a multidimensional array
% The |double| function can be used to convert a sumtensor to a multidimensional array.
% Similarly to the |full| expansion, this can use a prohibitive amount of
% memory for large-scale problems.
double(T)
%% Matricized Khatri-Rao product of a sumtensor
% The |mttkrp| function computes the Khatri-Rao product of a matricized tensor and a
% sumtensor. The required arguments are: a sumtensor X, a cell array of
% matrices U={U1,...,Um}, and a mode n. The cell array must consist of m matrices, 
% where m is the number of modes in X. The number of columns of these matrices
% should be constant, and number of rows of matrix Ui should match the dimension
% of the tensor X in mode i. The matricized Khatri-Rao product operation on 
% sumtensor distributes the operation to the summands of the sumtensor. 
% For details of this specific computation, see the mttkrp documentation 
% for a generic tensor. An example of the use of |mttkrp| follows.
U={eye(3), ones(3,3), randn(3,3)}; %<--The cell array of matrices
mttkrp(T,U,2)
%% Use innerprod to compute the inner product of a sumtensor
% The |innerprod| function computes the inner product of a sumtensor T and any type of 
% tensor. The operation is performed by distributing across each of the
% sumtensor's parts.
S = sptensor([1 1 1; 2 1 3; 3 2 2; 2 1 1], 1, [3,3,3]);
innerprod(T,S)
%% Use norm for compatibility with the other types of tensors.
% The |norm| function returns 0 and a warning when called on a sumtensor. 
% The procedure of computing the Frobenius norm of a sumtensor
% does not distribute across its parts, and hence is not supported for
% sumtensors. This default behavior is provided in order to ensure 
% compatibility of the sumtensor class with existing decomposition routines.
norm(T)
%%
% In order avoid this default behavior and return the Frobenius norm of a 
% sumtensor, it can be converted to a tensor using |full|.
norm(full(T))
%% Use CP-ALS to find a CP decomposition of a sumtensor
% One of the primary motivations for defining the |sumtensor| class is for
% efficient decomposition. In particular, when trying to find a CP
% decomposition of a tensor using alternating least squares, the
% subproblems can be efficiently created and solved using |mttkrp| and
% |innerprod|. Both of these operations can be performed more efficiently
% by exploiting extra structure in the tensors which form the sum, so the
% performance of |cp_als| is also improved. Consider the following example,
% where a |cp_als| is run on a sumtensor.
cp_als(T, 2)
%%
% It follows that in cases where $\mathcal{T}$ is too large for its full expansion to be 
% stored in memory, we may still be able find a CP decomposition by exploiting the 
% sumtensor structure. 
%%
% Note that the fit returned by cp_als is not correct for sumtensors, 
% because the norm operation is not supported.
%% Basic operations (plus) for sumtensors
% Sumtensors can be added to any other type of tensor. The result is a new
% sumtensor with the tensor appended to the parts of the original
% sumtensor. Note that the tensor is always appended, despite the order of
% the operation. 
T+S %<--S is appended to the parts of T
S+T %<--S is still the last part of T, despite order
%% Subscripted reference for sumtensors
% Subscripted reference can be used to return the individual parts of a 
% sumtensor.
T.part{1}
T.part{2}
