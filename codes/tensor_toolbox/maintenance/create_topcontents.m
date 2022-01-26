function create_topcontents
%CREATE_TOPCONTENTS Creates the top-level Contents.m file for toolbox

%% Open main contents file
fidmain = fopen('../Contents.m','w');
fprintf(fidmain,'%% Tensor Toolbox (Sandia National Labs)\n'); 
fprintf(fidmain,'%% Version 3.2.1 (R2021a) %s\n', date); 
fprintf(fidmain,'%% Tensor Toolbox for dense, sparse, and decomposed n-way arrays.\n'); 
fprintf(fidmain,'%% \n'); 
fprintf(fidmain,'%% Tensor Toolbox Classes:\n');
fprintf(fidmain,'%%   tensor     - Dense tensor.\n');
fprintf(fidmain,'%%   sptensor   - Sparse tensor.\n');
fprintf(fidmain,'%%   symtensor  - Symmetric tensor.\n');
fprintf(fidmain,'%%   ktensor    - Kruskal decomposed tensor.\n');
fprintf(fidmain,'%%   symktensor - Kruskal decomposed symmetric tensor.\n');
fprintf(fidmain,'%%   sumtensor  - Sum of different types of tensors.\n');
fprintf(fidmain,'%%   ttensor    - Tucker decomposed tensor.\n');
fprintf(fidmain,'%%   tenmat     - Tensor as matrix.\n');
fprintf(fidmain,'%%   sptenmat   - Sparse tensor as matrix.\n');
fprintf(fidmain,'%% \n'); 

%% Get contents of main directory
fprintf(fidmain,'%% Tensor Toolbox Functions:\n');
C = create_dircontents('..');
for i = 1:numel(C)
    fprintf(fidmain,'%%   %s\n',C{i});
end
fprintf(fidmain,'%%\n');
fprintf(fidmain,'%%   <a href="matlab:web(strcat(''file://'',fullfile(getfield(what(''tensor_toolbox''),''path''),''doc'',''html'',''index.html'')))">Documentation page for Tensor Toolbox</a>\n');
fprintf(fidmain,'%%\n');
fprintf(fidmain,'%%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>\n');



%% Close main contents file
fclose(fidmain);
