function disp(t, name)
%DISP Command window display for a symktensor.
%
%   DISP(T) displays a symmetric Kruskal tensor with no name.
%
%   DISP(T,NAME) display a symmetric Kruskal tensor with the given name.
%
%   See also DISP, SYMKTENSOR/DISPLAY, SYMKTENSOR
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if ~exist('name','var')
    name = 'ans';
end

fprintf('%s is a symktensor of order %d and dimension %d\n', name, t.m, size(t.u,1));
fprintf('\t%s.lambda = %s\n',name, ['[ ' num2str(t.lambda') ' ]'] );
fprintf('\t%s.U = \n', name);
output = tt_matrix2cellstr(t.u);
fprintf('\t\t%s\n',output{:});
