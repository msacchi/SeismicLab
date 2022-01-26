function t = subsasgn(t,s,b)
%SUBSASGN Subscripted assignment for tenmat.  
%
%   Examples 
%   X = tenmat(rand(3,4,2),1); 
%   X(1:2,1:2) = ones(2,2); <-- Calls SUBSASGN 
%
%   See also TENMAT, TENMAT/SUBSREF.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



switch s.type    
    case '()'
        [m n] = size(t.data);
	t.data(s.subs{:}) = b;
        if ~isequal([m n],size(t.data))
            error('Ambiguous change in size')
        end
    otherwise
        error('Invalid assignment for tenmat.')
end


