function disp(t,name)
%DISP Command window display of a sumtensor.
%
%   DISP(T) displays a sumtensor with no name.
%
%   DISP(T,NAME) display a sumtensor with the given name.
%
%   See also SUMTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if ~exist('name','var')
    name = 'ans';
end

if isempty(t.part)
    fprintf(1,'%s is an empty sumtensor\n', name);
    return;
end

fprintf(1,'%s is a sumtensor of size %s with %d parts\n', name, tt_size2str(size(t)), length(t.part));
for i = 1:length(t.part)
    subname = sprintf('%s.part{%d}',name,i);
    disp(t.part{i},subname);
end


