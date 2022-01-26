function s = tt_intvec2str(v)
%TT_INTVEC2STR Print integer vector to a string with brackets.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if isempty(v)
    s = sprintf('[]');
    return;
end

s = ['[ ' sprintf('%d ',v(1:end)) ']'];
