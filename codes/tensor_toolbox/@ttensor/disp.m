function disp(t,name)
%DISP Command window display of a ttensor.
%
%   DISP(T) displays a ttensor with no name.
%
%   DISP(T,NAME) display a ttensor with the given name.
%
%   See also TTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if ~exist('name','var')
    name = 'ans';
end

fprintf(1,'%s is a ttensor of size %s\n', name, tt_size2str(size(t)));
disp(t.core, sprintf('\t%s.core',name));

for j = 1 : ndims(t)
    fprintf('\t%s.U{%d} = \n', name, j);
    output = tt_matrix2cellstr(t.u{j});
    fprintf('\t\t%s\n',output{:});
end


