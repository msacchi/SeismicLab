function s = tt_size2str(sz)
%TT_SIZE2STR Convert size to a string that can be printed.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if isempty(sz)
    s = sprintf('[empty tensor]');
    return;
end

if numel(sz) == 1
    s = sprintf('%d',sz);
else
    s = [sprintf('%d x ',sz(1:end-1)) sprintf('%d', sz(end)) ];
end

