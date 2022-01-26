function disp(X, name) 
%DISP Command window display for a symtensor.
%
%   DISP(X) displays a symtensor with no name.
%
%   DISP(X,NAME) displays a tensor with the given name.
%
%   See also SYMTENSOR, SYMTENSOR/DISPLAY.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if ~exist('name','var')
    name = 'ans';
end

% preallocate
n=X.n;
m=X.m;
sz=length(X.val);
if sz==0   %empty array
   fprintf('%s is an empty symmetric tensor\n', name);
   return
end
output = cell(sz,1);

fprintf('%s is a symmetric tensor with %s modes of dimension %s\n',...
        name, num2str(X.m), num2str(X.n));

I=indices(X);

spc = floor(log10(max(double(I),[],1)))+1;
if numel(spc) == 1
    fmt = ['\t(%' num2str(spc(1)) 'd)%s'];
else
    fmt = ['\t(%' num2str(spc(1)) 'd,'];
    for i = 2:numel(spc)-1
        fmt = [fmt '%' num2str(spc(i)) 'd,'];
    end
    fmt = [fmt '%' num2str(spc(end)) 'd)%s'];
end
%%
% Get values out so that they look nice
savefmt = get(0,'FormatSpacing');
format compact
S = evalc('disp(X.val)');
set(0,'FormatSpacing',savefmt)
S = textscan(S,'%s','delimiter','\n','whitespace','');
S = S{1};
if ~isempty(strfind(S{1},'*'))
    fprintf('%s\n',S{1});
    S = S(2:end);
end
%%
for i = 1:sz
    output{i} = sprintf(fmt, I(i,:) ,S{i});
end
fprintf('%s\n',output{:});

end
    

