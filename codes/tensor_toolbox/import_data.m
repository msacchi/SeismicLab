function A = import_data(fname)
%IMPORT_DATA Import tensor-related data to a file.
%
%   A = IMPORT_DATA(FNAME) imports an object A from the file named FNAME.
%   The supported data types and formatting of the file are explained in
%   EXPORT_DATA. 
%
%   See also TENSOR, EXPORT_DATA
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%% Open file
fid = fopen(fname,'r');
if (fid == -1)
    error('Cannot open file %s',fname);
end

%% Get the type of object
type = import_type(fid);

%% Import the object

if strcmpi(type,'tensor')
    
    sz = import_size(fid);
    data = import_array(fid, prod(sz));   
    A = tensor(data, sz);
 
elseif strcmpi(type,'sptensor')
    
    sz = import_size(fid);
    nz = import_nnz(fid);
    [subs, vals] = import_sparse_array(fid, length(sz), nz);   
    A = sptensor(subs, vals, sz);
 
elseif strcmpi(type,'matrix') || strcmpi(type,'matrix')         

    sz = import_size(fid);
    data = import_array(fid, prod(sz));   
    A = reshape(data, sz);
    
elseif strcmpi(type,'ktensor')        

    sz = import_size(fid);
    r = import_rank(fid);
    lambda = import_array(fid, r);
    U = {};
    for n = 1:length(sz)
        line = fgets(fid);
        fac_type = import_type(fid);
        fac_sz = import_size(fid);
        fac_data = import_array(fid, prod(fac_sz));
        % row wise reshape
        fac = reshape(fac_data, fliplr(fac_sz))';
        U{n} = fac;
    end
    A = ktensor(lambda,U);
    
else   
    
    error('Invalid data type for export');    
    
end


%% Close file
fclose(fid);

function type = import_type(fid)
% Import IO data type
line = fgets(fid);
typelist = regexp(line, '\s+', 'split');
type = typelist(1);

function sz = import_size(fid)
% Import the size of something from a file
line = fgets(fid);
n = sscanf(line, '%d');
line = fgets(fid);
sz = sscanf(line, '%d');
sz = sz';
if (size(sz,2) ~= n)
    error('Imported dimensions are not of expected size');
end

function nz = import_nnz(fid)
% Import the size of something from a file
line = fgets(fid);
nz = sscanf(line, '%d');

function r = import_rank(fid)
% Import the rank of something from a file
line = fgets(fid);
r = sscanf(line, '%d');

function data = import_array(fid, n)
% Import dense data that supports numel and linear indexing
data = fscanf(fid, '%e', n);

function [subs, vals] = import_sparse_array(fid, n, nz)
% Import sparse data subs and vals from coordinate format data
data = textscan(fid,[repmat('%f',1,n) '%n']);
subs = cell2mat(data(1:n));
vals = data{n+1};
if (size(subs,1) ~= nz)
    error('Imported nonzeros are not of expected size');
end
