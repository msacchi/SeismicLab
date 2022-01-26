%% Symmetric Tensors
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="tensor_types.html">Tensor Types</a> 
% &#62;&#62; <a href="symtensor_doc.html">Symmetric Tensors</a>
% </p>
% </html>
%
% A symmetric tensor is a tensor that is invariant under all permutations
% of its modes.  Because many of the elements of a symmetric tensor are
% guaranteed to be equal, we can simplify the storage requirements by only
% storing the unique values of the symmetric tensor. There are
% ${n+m-1}\choose{m}$ such values for an m-way tensor of dimension n. 
% The |symtensor| class is designed to take advantage of this symmetric
% structure by only storing the unique values of the tensor. 

%% Definition of a symmetric tensor
% A symmetric tensor is invariant under any permutation of the indices.
% Here is a small example. The |issymmetric| function checks symmetry of a
% dense tensor.
T(:,:,1) = [1 2; 2 3]; T(:,:,2)= [2 3; 3 6]; 
T = tensor(T) 
issymmetric(T)

%% Creating a symtensor from a symmetric tensor
% We can construct a |symtensor| object from a symmetric tensor. This
% object only stores the unique entries of the tensor. For the 2 x 2 x 2
% tensor, this means there are only four unique entries. Everything else
% comes from permuting the indices of those four entries.
S = symtensor(T) 

%% Unique entries of a tensor
% Note from TGK: This needs to be added. It should have some discussion of
% all the return values from indices. What is a monomial description, etc.
[I,C,W,Q] = indices(S)

%% Creating a symtensor from a nonsymmetric tensor
% A symmetric tensors can be created from the symmetrization of
% nonsymmetric tensor so long as it is the same size in every mode. 
% If the input is not symmetric, it is symmetrized by creating an average
% of elements in the same permutation class. For instance, this example
% starts with a nonsymmetric tensor and symmetrizes it:
T2 = tensor([1:8],[2 2 2])
S2 = symtensor(T2)

%% 
% Converting the symtensor back to a generic tensor is equivalent to
% running |symmetrize| on the original tensor. In the following example,
% the full command converts a symtensor to a tensor.
full(S2)
symmetrize(T2)

%% Create an all ones symtensor
% The first argument is the generating function, the second argument is the
% number of modes, and the third argument is the size of each mode.
S3 = symtensor(@ones, 3, 2)

%% Create a random symtensor
S4 = symtensor(@randn, 3, 2)

%% Using a generating function to populate a symmetric tensor
% In general, a symmetric tensor can also have its entries created by any
% generating function. This is done by passing a function handle, the
% number of modes, and the dimension. The function is expected to take a
% two-dimension size as input and return a matrix of that shape. In fact,
% the second argument to the function will always be 1. 

% For example, we can also declare a binary symmetric tensor as follows:
S5 = symtensor(@(x,y) double(rand(x,y)>.25), 3, 3)

%% Use ndims and size to get the size of a symmetric tensor
ndims(S) %<-- Number of modes of the symmetric tensor

%%
size(S) %<-- Size of a symmetric tensor

%% Use full to convert a symmetric tensor to a multidimensional array
full(S) %<-- Converts from a symmetric tensor to a tensor

%% Subscripted reference of a symmetric tensor
% Subindex notation extracts the tensor value. 
S(1,2,2)
S(2,1,2) %<-- Equal to above, by symmetry

%%
% This works the same as applying it to the full tensor.
T(1,2,2)
T(2,1,2)

%%
% Multiple indices can be queried by combining these indices into the rows 
% of a matrix. Consider the following example, which returns a vector 
% consisting of the values of S at indices indicated by the rows of the
% input matrix.
S([1 2 1; 2 1 2])

%%
% Single indices are interpretted as an index into the unique value array,
% which is stored with respect to increasing indices. This is very
% different than using linear indexing on the full tensor.
S(3) %<- Third unique entry corresponding to (1,2,2)
S(4) %<- Fourth unique entry, corresponding to (2,2,2)
T(3) %<- Third entry in the tensor, i.e., (1,2,1) = (1,1,2)
T(4) %<- Fourth entry in the tensor, i.e., (2,2,1) = (1,2,2)

%%
% Mulitple entries can be obtained at once as well.
S([3:4]')

%% 
%% Subscripted assignment
% Symmetric tensors also support subscripted assignment. Either linear or
% subindex notation is valid. Multiple values can be assigned the same
% quantity, but assigning a subset of a symmetric tensor from a 
% multidimensional arrays, tensor, or symtensor is not allowed.
S5(1) = 7 %<-- Linear indexing
S5(2,1,2) = 6 %<-- Subindex indexing
%%
% Symmetric tensors do not support enlargement with the assignment
% operator, so assigning a value to an index other than those which have
% already been declared produces an error.
%% Basic operations (plus, minus, and, or, etc.) on a symmetric tensor
% The tensor object supports many basic operations, illustrated here.
A = symtensor(@(x,y) rand(x,y)>.5, 3, 2)
B = symtensor(@(x,y) rand(x,y)>.5, 3, 2)
%%
A==B %<-- Calls eq.
%%
A<B %<-- Calls lt.
%%
A.*B %<-- Calls times. (elementwise multiply)
%%
5*A %<-- Calls mtimes. (scalar multiply)
%%
% The symtensor class supports the following additional MATLABbinary
% operations: and, or, xor, neq, gt, ge, le, plus, minus, power, ldivide,
% and rdivide. Supported unary operations include: not, uplus, uminus.

%% Using |symtenfun| for elementwise operations on one or more symmetric tensors
% The function |symtenfun| applies a function to a number of symmetric
% symtensors. This function mirrors the capability of |tenfun| for tensors.
tenfun(@min, S, S2, S+1) %<-- Symtensor formed from elementwise minimization
%%
tenfun(@(x)(x-1.5),S) %<-- Subtract 1.5 from each element of B
