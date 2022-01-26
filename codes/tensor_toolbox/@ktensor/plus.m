function C = plus(A,B)
%PLUS Binary addition for ktensor.
%
%   C = PLUS(A,B) adds two ktensors of the same size, and the
%   result is a ktensor of the same size.
%
%   See also KTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if isa(B,'sumtensor') %If the 2nd argument is a sumtensor
    C = plus(B,A); %Call plus for sumtensor
    return
end

if (isa(A,'ktensor') && isa(B,'ktensor'))    

    if ~isequal(size(A),size(B))
	error('Tensor size mismatch.')
    end

    lambda = [A.lambda; B.lambda];    
    M = ndims(A);
    u = cell(M,1);
    for m = 1 : M
        u{m} = [A.u{m} B.u{m}];
    end 
    C = ktensor(lambda, u);
    return;
end

error('Use plus(full(A),full(B)).');
