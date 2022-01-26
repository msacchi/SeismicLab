% Tests for symtensor class
classdef Test_Symtensor < matlab.unittest.TestCase
    
    properties (TestParameter)
        m = struct( 'three', 3, 'one', 1 );
        n = struct( 'ten', 10, 'three', 3, 'one', 1);
        g = struct( 'rand', @rand, 'zeros', @zeros, 'ones', @ones);
        
        bf = struct('and', @and, 'eq', @eq, 'ge', @ge, 'gt', @gt, 'ldivide', @ldivide, 'le', @le, 'lt', @lt, 'minus', @minus, 'ne', @ne, 'or', @or, 'plus', @plus, 'power', @power, 'rdivide', @rdivide, 'times', @times, 'xor', @xor);
        uf = struct('not', @not, 'uminus', @uminus, 'uplus', @uplus);
        sfr = struct('mtimes', @mtimes, 'mrdivide', @mrdivide);
        sfl = struct('mldivide', @mldivide, 'mtimes', @mtimes);
    end
    
    methods (Test)
        
        function ConstructBySym(testCase,m,n,g)
            % Test construction by symmetrization
            % S = SYMTENSOR(X)
            % [S,I] = SYMTENSOR(X)
            lsz = n * ones(1,m);
            T = tensor(feval(g,[lsz 1]),lsz);
            [X,I] = symtensor(T);
            testCase.verifyClass(X, 'symtensor');
            testCase.verifyEqual(X.m, m);
            testCase.verifyEqual(X.n, n);
            testCase.verifyEqual(full(X), symmetrize(T));
            testCase.verifyEqual(size(I), [nchoosek(m+n-1,m) m]);
            % Also test indices function while we're here...
            I2 = indices(X);
            testCase.verifyEqual(I, I2);
            % Test some other stuff too
            testCase.verifyTrue(issymmetric(X)); % Also testing issymmetric
            testCase.verifyEqual(ndims(X), m);
            testCase.verifyEqual(size(X), n*ones(1,m));
            testCase.verifyEqual(size(X,randi(m)), n);
            % Test isequal too
            Y = symtensor(X.vals, m, n);
            testCase.verifyTrue(isequal(X,Y));
            Y(1) = Y(1) - 1;
            testCase.verifyFalse(isequal(X,Y));
        end
        
        function ConstructByCopy(testCase)
            % Test copy constructor
            % S = SYMTENSOR(S0)
            X = symtensor(@rand,4,3);
            Y = X;
            testCase.verifyClass(Y, 'symtensor');
            testCase.verifyEqual(X,Y);
        end
        
        function ConstructByVal(testCase, m, n, g)
            % Test construct from distinct values
            %   S = SYMTENSOR(VALS,M,N)
            vsz = nchoosek(m+n-1,m);
            vals = feval(g,[vsz 1]);
            X = symtensor(vals,m,n);
            Y = symtensor(vals',m,n);
            testCase.verifyClass(X, 'symtensor');
            testCase.verifyClass(Y, 'symtensor');
            testCase.verifyEqual(X,Y);
            % Test some other stuff too
            testCase.verifyTrue(issymmetric(X)); % Also testing issymmetric
            testCase.verifyEqual(ndims(X), m);
            testCase.verifyEqual(size(X), n*ones(1,m));
            testCase.verifyEqual(size(X,randi(m)), n);
            % Test isequal too
            Y = symtensor(X.vals, m, n);
            testCase.verifyTrue(isequal(X,Y));
            Y(1) = Y(1) - 1;
            testCase.verifyFalse(isequal(X,Y));
        end
        
        function ConstructByGen(testCase, m, n, g)
            X = symtensor(g, m, n);
            testCase.verifyClass(X, 'symtensor');
            testCase.verifyEqual(X.m, m);
            testCase.verifyEqual(X.n, n);
            % Test some other stuff too
            testCase.verifyTrue(issymmetric(X)); % Also testing issymmetric
            testCase.verifyEqual(ndims(X), m);
            testCase.verifyEqual(size(X), n*ones(1,m));
            testCase.verifyEqual(size(X,randi(m)), n);
            % Test isequal too
            Y = symtensor(X.vals, m, n);
            testCase.verifyTrue(isequal(X,Y));
            Y(1) = Y(1) - 1;
            testCase.verifyFalse(isequal(X,Y));
        end
        
        function BinaryFuncs(testCase, bf)
            mlocal = 4;
            nlocal = 3;
            X = symtensor(@ones,mlocal,nlocal); % 15 distinct elements
            Y = X+1;
            Y((1:5)') = 0;
            Z1 = feval(bf, X, Y);
            Z2 = symtensor(feval(bf, X.vals, Y.vals), mlocal, nlocal);
            testCase.verifyEqual(Z1,Z2);
        end
        
        function BinaryFuncsWithScalar(testCase, bf)
            mlocal = 4;
            nlocal = 3;
            X = symtensor(@ones, mlocal, nlocal); % 15 distinct elements
            Z1 = feval(bf, X, 5);
            Z2 = symtensor(feval(bf, X.vals, 5), mlocal, nlocal);
            testCase.verifyEqual(Z1,Z2);
        end
        
        function UnaryFuncs(testCase,uf)
            mlocal = 3;
            nlocal = 5;
            X = symtensor(@rand, mlocal, nlocal);
            X((1:5)') = X((1:5)') > .5;
            Z1 = feval(uf,X);
            Z2 = symtensor(feval(uf,X.vals), mlocal, nlocal);
            testCase.verifyEqual(Z1,Z2);
        end
        
        function ScalarFuncsRight(testCase,sfr)
            mlocal = 4;
            nlocal = 3;
            X = symtensor(@rand, mlocal, nlocal); % 15 distinct elements
            Z1 = feval(sfr, X, 5);
            Z2 = symtensor(feval(sfr, X.vals, 5), mlocal, nlocal);
            testCase.verifyEqual(Z1,Z2);
        end
        
        function ScalarFuncsLeft(testCase,sfl)
            mlocal = 4;
            nlocal = 3;
            X = symtensor(@rand, mlocal, nlocal); % 15 distinct elements
            Z1 = feval(sfl, 5, X);
            Z2 = symtensor(feval(sfl, 5, X.vals), mlocal, nlocal);
            testCase.verifyEqual(Z1,Z2);
        end
        
        function CheckFull(testCase)
            X=tensor(rand([3,3,3]));
            Y=symtensor(X);
            F=full(Y);
            testCase.verifyEqual(symmetrize(X),F);
            Y=symtensor(@ones,3,3);
            F=full(Y);
            testCase.verifyEqual(F,tensor(ones([3,3,3])));
        end
        
        function SubsRef(testCase, m, n)
            p = nchoosek(m+n-1,m);
            X = symtensor(1:p, m, n);
            testCase.verifyEqual(X.val, (1:p)');
            testCase.verifyEqual(X.m, m);
            testCase.verifyEqual(X.n, n);
            %% Linear indexing into val array
            q = randi(p);
            % Single
            testCase.verifyEqual(X(q),q);
            % Range
            testCase.verifyEqual(X((1:q)'),(1:q)');
            r = (min(2,p):min(4,p))';
            testCase.verifyEqual(X(r),r);
            % List
            q = randi(p,25,1);
            testCase.verifyEqual(X(q),q);
            %% Subscripts
            % Single
            s = randi(n,1,m);
            ssrt = sort(s,2);
            testCase.verifyEqual(X(s),X(ssrt));
            xsubs = indices(X);
            [~,locx] = ismember(ssrt,xsubs,'rows');
            testCase.verifyEqual(X(s),X(locx));
            % List
            s = randi(n,5,m);
            ssrt = sort(s,2);
            testCase.verifyEqual(X(s),X(ssrt));
            xsubs = indices(X);
            [~,locx] = ismember(ssrt,xsubs,'rows');
            testCase.verifyEqual(X(s),X(locx));
        end
        
        function SubsAsgn(testCase, m, n)
            p = nchoosek(m+n-1,m);
            X = symtensor(1:p, m, n);
            % Assignment of all values
            X(:) = (p:-1:1)';
            testCase.verifyEqual(X.val, (p:-1:1)');
            X.val = (1:p)';
            testCase.verifyEqual(X.val, (1:p)');
            % Assignment of single entry using linear index
            idx = randi(p);
            X(idx) = -X(idx);
            newvals = (1:p)';
            newvals(idx) = -idx;
            testCase.verifyEqual(X.val, newvals);
            % Assignment of multiple entries using linear indices
            idx = unique(randi(p,10,1));
            newvals = -1*(1:length(idx))';
            X(idx) = newvals;
            testCase.verifyEqual(X(idx),newvals);
            % Assignment of single entry using subscripts
            s = randi(n,1,m);
            X(s) = 17;
            testCase.verifyEqual(X(s),17);
            % Assignment of multiple entries using subscripts
            % (Need to be careful here because we can have repeat entries in s
            % that are not obvious until sorted due to different subscripts for
            % the same unique element.)
            s = randi(n,5,m);
            ssrt = sort(s,2);
            [ssrt,locs] = unique(ssrt,'rows');
            s = s(locs,:);
            newvals = 10*(1:size(s,1))';
            X(s) = newvals;
            testCase.verifyEqual(X(ssrt),newvals);
        end
        
        function TenFun(testCase, m, n, g)
            X = symtensor(g, m, n);
            Xf = full(X);
            fh = @(x) x + 1;
            Z  = tenfun(fh, X );
            Zf = tenfun(fh, Xf);
            testCase.verifyEqual(Z, symtensor(Zf));
            
            fh = @eq;
            Z  = tenfun(fh, X , 1);
            Zf = tenfun(fh, Xf, 1);
            testCase.verifyEqual(Z, symtensor(Zf));
            
            Y = symtensor(g, m, n);
            Yf = full(Y);
            fh = @plus;
            Z  = tenfun(fh, X , Y);
            Zf = tenfun(fh, Xf, Yf);
            testCase.verifyEqual(Z, symtensor(Zf));
            
            W = symtensor(g, m, n);
            Wf = full(W);
            fh = @max;
            Z  = tenfun(fh, W, X, Y);
            Zf = tenfun(fh, Wf, Xf, Yf);
            testCase.verifyEqual(Z, symtensor(Zf));
        end
        
    end
end


