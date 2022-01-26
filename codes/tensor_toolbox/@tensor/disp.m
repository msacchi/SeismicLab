function disp(X,name)
%DISP Command window display of a tensor.
%
%   DISP(X) displays a tensor with no name.
%
%   DISP(X,NAME) displays a tensor with the given name.
%
%   See also TENSOR, TENSOR/DISPLAY.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if ~exist('name','var')
    name = 'ans';
end

fprintf(1,'%s is a tensor of size %s\n',name,tt_size2str(X.size));

if isempty(X.data)
    fprintf(1,'\t%s = []\n',name);
    return
end

s = shiftdim(num2cell(X.data,1:2),2);

for i = 1:numel(s)
    fprintf('\t%s',name);
    if ndims(X) == 1
        fprintf('(:)');
    elseif ndims(X) == 2
        fprintf('(:,:)');
    elseif ndims(X) > 2
        fprintf('(:,:');
        fprintf(',%d',tt_ind2sub(X.size(3:end),i));
        fprintf(')');
    end
    fprintf(' = \n');
    output = tt_matrix2cellstr(s{i});
    fprintf('\t%s\n',output{:});
end

end
