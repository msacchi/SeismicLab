function C = minus(A,B)
%MINUS Binary subtraction for ktensor.  
%
%   C = MINUS(A,B) computes C = A - B.  A and B must both be ktensors
%   and have the same size, and the result is another ktensor of the
%   same size.
%
%   C = MINUS(A,B) is called for the syntax 'A - B' when A or B is a
%   ktensor.
%
%   See also KTENSOR, SIZE, ISEQUAL.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if (isa(A,'ktensor') && isa(B,'ktensor'))    

    if ~isequal(size(A),size(B))
	error('Tensor size mismatch.')
    end

    lambda = [A.lambda; -B.lambda];    
    M = ndims(A);
    u = cell(M,1);
    for m = 1 : M
        u{m} = [A.u{m} B.u{m}];
    end 
    C = ktensor(lambda,u);
    return;

end

error('Use minus(full(A),full(B)).');



