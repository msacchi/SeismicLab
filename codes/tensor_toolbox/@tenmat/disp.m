function disp(t,name)
%DISP Command window display of a matricized tensor (tenmat).
%
%   DISP(T) displays a tensor as matrix with no name.
%
%   DISP(T,NAME) display a tensor as matrix with the given name.
%
%   See also TENMAT, TENMAT/DISPLAY.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if ~exist('name','var')
    name = 'ans';
end
    
fprintf('%s is a matrix corresponding to a tensor of size %s\n',...
        name,tt_size2str(t.tsize));
fprintf('\t%s.rindices = %s (modes of tensor corresponding to rows)\n',...
        name,['[ ' num2str(t.rindices) ' ]']);
fprintf('\t%s.cindices = %s (modes of tensor corresponding to columns)\n',...
        name,['[ ' num2str(t.cindices) ' ]']);

if isempty(t.data)
    fprintf('\t%s.data = []\n',name);
else
    fprintf('\t%s.data = \n',name);
    output = tt_matrix2cellstr(t.data);
    fprintf('\t\t%s\n',output{:});
end    


