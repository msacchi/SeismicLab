function a = reshape(a,new_size,old_modes)
%RESHAPE Reshape sparse tensor.
%   
%   RESHAPE(X,SIZ) reshapes the sparse tensor to the given size. PROD(SIZ)
%   must be the same as PROD(SIZE(X)).
%
%   RESHAPE(X,SIZ,MODES) reshapes only the specifies modes and appends the
%   new reshaped modes to the end of the indices.
%
%   See also SPTENSOR, SPTENSOR/PERMUTE, RESHAPE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if ~exist('old_modes','var')
    old_modes = 1:ndims(a);
    keep_modes = [];
else
    keep_modes = setdiff(1:ndims(a),old_modes);
end
old_size = a.size(old_modes);
keep_size = a.size(keep_modes);


if prod(new_size) ~= prod(old_size)
    error('prod(SIZ) must be the same size of prod(SIZE(X,MODES))');
end

if isempty(a.subs)
    a.size = [keep_size new_size];
    a.subs = [];
else
    inds = tt_sub2ind(old_size,a.subs(:,old_modes));
    new_subs = tt_ind2sub(new_size,inds);

    a.size = [keep_size new_size];
    a.subs = [a.subs(:,keep_modes) new_subs];
end

