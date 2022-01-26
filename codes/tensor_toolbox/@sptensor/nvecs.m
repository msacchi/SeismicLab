function u = nvecs(t,n,r,opts)
%NVECS Compute the leading mode-n vectors for a sparse tensor.
%
%   U = NVECS(X,n,r) computes the r leading eigenvalues of Xn*Xn'
%   (where Xn is the mode-n matricization of X), which provides
%   information about the mode-n fibers. In two-dimensions, the r
%   leading mode-1 vectors are the same as the r left singular vectors
%   and the r leading mode-2 vectors are the same as the r right
%   singular vectors.
%
%   U = NVECS(X,n,r,OPTS) specifies options:
%   OPTS.eigsopts: options passed to the EIGS routine [struct('disp',0)]
%   OPTS.flipsign: make each column's largest element positive [true]
%
%   Examples
%   S = sptensor([3 3 3; 1 3 2; 1 1 3], 1, [3,3,3]); %<--Declare an sptensor
%   nvecs(S,3,2)
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','nvecs_doc.html')))">Documentation page for n-vecs</a>
%
%   See also SPTENSOR, SPTENMAT, EIGS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if ~exist('opts','var')
    opts = struct;
end


if isfield(opts,'eigsopts')
    eigsopts = opts.eigsopts;
else
    eigsopts.disp = 0;
end

tnt = double(sptenmat(t,n,'t'));
y = tnt' * tnt;
opts.disp = 0;
[u,d] = eigs(y,r,'LM',eigsopts);

%tn = sptenmat(t,n);
%[u,d] = eigs(@(x)aatx(tn,x), size(t,n), r, 'LM', eigsopts);

if isfield(opts,'flipsign') 
    flipsign = opts.flipsign;
else
    flipsign = true;
end
    
if flipsign
    % Make the largest magnitude element be positive
    [val,loc] = max(abs(u));
    for i = 1:r
        if u(loc(i),i) < 0
            u(:,i) = u(:,i) * -1;
        end
    end
end
