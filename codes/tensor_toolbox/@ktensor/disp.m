function disp(t, name)
%DISP Command window display for a ktensor.
%
%   DISP(T) displays a Kruskal tensor with no name.
%
%   DISP(T,NAME) display a Kruskal tensor with the given name.
%
%   See also DISP, KTENSOR/DISPLAY, KTENSOR
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if ~exist('name','var')
    name = 'ans';
end

fprintf('%s is a ktensor of size %s\n', name, tt_size2str(size(t)));
output = tt_matrix2cellstr(t.lambda');
fprintf('\t%s.lambda = \n',name);
fprintf('\t\t%s\n',output{:});

if (ndims(t) > 0)
    for j = 1 : ndims(t)
        fprintf('\t%s.U{%d} = \n', name, j);
        output = tt_matrix2cellstr(t.u{j});
        fprintf('\t\t%s\n',output{:});
    end
end
