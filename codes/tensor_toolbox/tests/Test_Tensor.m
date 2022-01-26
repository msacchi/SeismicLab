% Tests for tensor class
classdef Test_Tensor < matlab.unittest.TestCase
    
    properties (TestParameter)
        nd = struct( 'three', 3, 'one', 1 );
        maxdim = struct( 'ten', 10, 'three', 3, 'one', 1 );
        gen = struct( 'rand', @rand, 'zeros', @zeros, 'ones', @ones, 'randn', @randn );
        nx = struct( 'one', 1, 'ten', 10, 'hundred', 100 );
        szs = struct( 'fourthreetwo', [4 3 2], 'threetwoone', [3 2 1], 'threetwo', [3 2], 'three', 3, 'one', 1, 'emtpy', []);
        
        bf = struct('and', @and, 'eq', @eq, 'ge', @ge, 'gt', @gt, 'ldivide', @ldivide, 'le', @le, 'lt', @lt, 'minus', @minus, 'ne', @ne, 'or', @or, 'plus', @plus, 'power', @power, 'rdivide', @rdivide, 'times', @times, 'xor', @xor);
        uf = struct('not', @not, 'uminus', @uminus, 'uplus', @uplus);
        sfr = struct('mtimes', @mtimes, 'mrdivide', @mrdivide);
        sfl = struct('mldivide', @mldivide, 'mtimes', @mtimes);
        
        symver = struct('new', 0, 'old', 1);
        issymver = struct('new', 0', 'old', 1);
        
    end
    
    methods (Test)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- CONSTRUCTOR TESTS ---
        % Testing the constructor using various generators and also the
        % (implicit) copy constructor. Also including some other tests that
        % are simple: double, ndims, full, isequal, norm, nnz, size.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        function Construct(testCase, nd, maxdim, gen)
            sz = randi(maxdim, 1, nd);
            % --- Case I: Construct from (Vector) Data ---
            s = rng; % Save state of random number generator
            Xdata = gen([sz 1]);
            X = tensor(Xdata, sz);
            % --- Case II: Use Generator Function ---
            rng(s); % Reset random number generator to same state as above
            Y = tensor(gen,sz); % Should be the same as X
            % --- Case III: Copy ---
            Z = X;
            % --- Checks ---
            testCase.verifyClass(X, 'tensor');
            testCase.verifyClass(Y, 'tensor');
            testCase.verifyClass(Z, 'tensor');
            testCase.verifyEqual(X, Y);
            testCase.verifyEqual(X, Z);
            testCase.verifyEqual(ndims(X), nd);
            testCase.verifyEqual(size(X), sz);
            for i = 1:nd;
                testCase.verifyEqual(size(X,i), sz(i));
            end
            % --- More Tests ---
            testCase.verifyEqual(double(X), Xdata); % double
            testCase.verifyEqual(norm(X), sqrt(sum(Xdata(:).^2)),'RelTol',1e-15); % norm
            testCase.verifyEqual(X, full(X)); % full
            testCase.verifyTrue(isequal(X,Y)); % isequal
            testCase.verifyEqual(nnz(X), nnz(X.data(:)));
        end
        
        function ConstructAlt(testCase,szs)
            X = tenrand(szs);
            testCase.verifyEqual(size(X),szs);
        end
        
        function ConstructBadSize(testCase)
            testCase.verifyError(@()eval('X = tensor(1:12, [4 3 2]);'), ?MException);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- FUNCTION TESTS ---
        % Testing binary and unary functions to make sure they work as
        % expected. In particular, the binary functions should accept
        % scalar inputs. The matrix functions (mtimes, mrdivide, mldivide)
        % work with scalars only.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        function BinaryFuncs(testCase, bf)
            sz = [4 3 2 4];
            X = tensor(@ones, sz);
            Y = X+1;
            Y((1:5)') = 0;
            Z1 = feval(bf, X, Y);
            Z2 = tensor(feval(bf, X.data, Y.data), sz);
            testCase.verifyEqual(Z1,Z2);
        end
        
         
        function BinaryFuncErrors(testCase,bf) %check binary functions
            X = tensor(rand([4 3 2]));
            Y = tensor(rand([3 4 2]));
            % Should throw an exception if the sizes don't match
            testCase.verifyError(@()bf(X,Y),?MException);
            % Should throw an exception if only a single argument
            testCase.verifyError(@()bf(X), ?MException);
            % Should throw an exception for too many arguments
            testCase.verifyError(@()bf(X,X,X), ?MException);
        end

        function And(testCase, szs)
            X = tenrand(szs);
            Y = X.data;
            testCase.verifyEqual(double(X&X), double(Y~=0));
        end
        
        function UnaryFuncs(testCase,uf)
            sz = [4 3 2 4];
            X = tensor(@rand, sz);
            X((1:5)') = X((1:5)') > .5;
            Z1 = feval(uf,X);
            Z2 = tensor(feval(uf,X.data), sz);
            testCase.verifyEqual(Z1,Z2);
        end
        
        function ScalarFuncsRight(testCase,sfr)
            sz = [4 3 2 4];
            X = tensor(@rand, sz);
            Z1 = feval(sfr, X, 5);
            Z2 = tensor(feval(sfr, X.data, 5), sz);
            testCase.verifyEqual(Z1,Z2);
        end
        
        function ScalarFuncsLeft(testCase,sfl)
            sz = [7 4 1];
            X = tensor(@rand, sz);
            Z1 = feval(sfl, 5, X);
            Z2 = tensor(feval(sfl, 5, X.data), sz);
            testCase.verifyEqual(Z1,Z2);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- SUBSCRIPTED REFERENCE AND ASSIGNMENT ---
        % Testing the various ways of referencing a tensor, including dot,
        % subscripts, and linear indices. Tensor support passing an array
        % of linear indices or a matrix of subscripts with one per row.
        % Tensors also support extraction of a subtensor using ranges.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function SubsRefDot(testCase, nd ,maxdim, gen)
            sz = randi(maxdim, 1, nd);
            Xdata = gen([sz 1]);
            X = tensor(Xdata, sz);
            testCase.verifyEqual(X.data, Xdata);
            testCase.verifyEqual(X.size, sz);
        end
        
        function SubsRefEntry(testCase, nd, maxdim, gen)
            sz = randi(maxdim, 1, nd);
            idx = randi(prod(sz));
            sub = tt_ind2sub(sz, idx);
            
            szstr = strjoin(arrayfun(@num2str, sz, 'UniformOutput', false),',');
            substr = strjoin(arrayfun(@num2str, sub, 'UniformOutput', false),',');
            
            Xdata = gen([sz 1]);
            X = tensor(Xdata, sz);
            
            dgnstc = sprintf('Failure for size=%s, entry=%s', szstr, substr);
            testCase.verifyEqual(Xdata(idx), X(idx), dgnstc);
            testCase.verifyEqual( eval(['X(' substr ')']), X(idx), dgnstc);
        end
        
        function SubsRefColons(testCase, nd, maxdim, gen)
            sz = randi(maxdim, 1, nd);
            szstr = strjoin(arrayfun(@num2str, sz, 'UniformOutput', false),',');
            
            Xdata = gen([sz 1]);
            X = tensor(Xdata, sz);             %#ok<NASGU>
            
            iscolon = rand(nd,1) > 0.5;
            iscolon(randi(nd)) = 1;
            
            substrs = cell(nd,1);
            for k = 1:nd
                if iscolon(k)
                    substrs{k} = ':';
                else
                    substrs{k} = num2str(randi(sz(k)));
                end
            end
            substr = strjoin(substrs,',');
            
            Y1 = eval(['X(' substr ')']);
            Y2 = eval(['Xdata(' substr ')']);
            dgnstc = sprintf('Failure for size=%s, entry=%s', szstr, substr);
            testCase.verifyEqual(Y1.data(:), Y2(:), dgnstc);
            testCase.verifyEqual(size(Y1), sz(iscolon), dgnstc);
        end
        
        function SubsRefRanges(testCase, nd, maxdim, gen)
            sz = randi(maxdim, 1, nd);
            szstr = strjoin(arrayfun(@num2str, sz, 'UniformOutput', false),',');
            
            Xdata = gen([sz 1]);
            X = tensor(Xdata, sz);             %#ok<NASGU>
            
            iscolon = rand(nd,1) > 0.5;
            iscolon(randi(nd)) = 1;
            
            newsz = sz;
            substrs = cell(nd,1);
            for k = 1:nd
                if iscolon(k)
                    if sz(k) == 1
                        substrs{k} = ':';
                    else
                        lower = randi(sz(k)-1);
                        upper = lower+randi(sz(k)-lower);
                        substrs{k} = sprintf('%d:%d', lower, upper);
                        newsz(k) = upper - lower + 1;
                    end
                else
                    substrs{k} = num2str(randi(sz(k)));
                end
            end
            substr = strjoin(substrs,',');
            
            Y1 = eval(['X(' substr ')']);
            Y2 = eval(['Xdata(' substr ')']);
            dgnstc = sprintf('Failure for size=%s, entry=%s', szstr, substr);
            testCase.verifyEqual(Y1.data(:), Y2(:), dgnstc);
            testCase.verifyEqual(size(Y1), newsz(iscolon), dgnstc);
        end
        
        function SubsRefList(testCase, nd, maxdim, nx)
            sz = randi(maxdim, 1, nd);
            X = tensor(@rand, sz);
            idx = randi(prod(sz),nx,1); % List of elements to extract
            sub = tt_ind2sub(sz, idx);
            L0 = squeeze(X.data(idx));
            if ~iscolumn(L0)
                L0 = L0';
            end
            L1 = X(idx, 'extract');
            testCase.verifyEqual(L0,L1);
            L2 = X(sub, 'extract');
            testCase.verifyEqual(L0,L2);
            if nd > 1
                L3 = X(idx);
                testCase.verifyEqual(L0,L3);
                L4 = X(sub);
                testCase.verifyEqual(L0,L4);
            end
        end
        
        
        function SubsAsgnElement(testCase, nd, maxdim)
            sz = randi(maxdim, 1, nd);
            szstr = strjoin(arrayfun(@num2str, sz, 'UniformOutput', false),',');
            idx = randi(prod(sz));
            sub = tt_ind2sub(sz, idx);
            substr = strjoin(arrayfun(@num2str, sub, 'UniformOutput', false),',');
            dstr = sprintf('Failure sz=%s, idx=%d, sub=%s\n', szstr, idx, substr);
            
            % Using linear index
            X = tensor(@zeros, sz);
            X(idx) = 1;
            testCase.verifyEqual(X(idx), 1, dstr);
            if prod(sz) > 1
                rng = [1:idx-1, idx+1:prod(sz)]';
                testCase.verifyEqual(X(rng,'extract'), zeros(prod(sz)-1,1),dstr);
            end
            
            % Repeat using subscript
            X = tensor(@zeros, sz);
            estr = sprintf('X( %s ) = 1;', substr);
            eval(estr);
            testCase.verifyEqual(X(idx), 1, dstr);
            testCase.verifyEqual(nnz(X), 1, dstr);
        end
        
        function SubsAsgnGrowSize(testCase, nd, maxdim)
            sub = randi(maxdim, 1, nd);
            substr = strjoin(arrayfun(@num2str, sub, 'UniformOutput', false),',');
            X = tensor;
            estr = sprintf('X( %s ) = 1;', substr);
            eval(estr);
            testCase.verifyEqual(size(X), sub);
            testCase.verifyEqual(X.data(end), 1);
            testCase.verifyEqual(nnz(X), 1);
        end
        
        function SubsAsgnGrowOrder(testCase, nd, maxdim)
            sz = randi(maxdim, 1, nd);
            sub = [tt_ind2sub(sz, randi(prod(sz))) 1];
            newsz = [sz 1];
            idx = sub2ind(newsz, sub);
            substr = strjoin(arrayfun(@num2str, sub, 'UniformOutput', false),',');
            X = tensor(@zeros,sz);
            estr = sprintf('X( %s ) = 1;', substr);
            eval(estr);
            testCase.verifyEqual(size(X), newsz);
            testCase.verifyEqual(X(idx), 1);
            testCase.verifyEqual(ndims(X), nd+1);
            testCase.verifyEqual(nnz(X), 1);
        end
        
        function SubsAsgnListValue(testCase, nd, maxdim, nx)
            sz = randi(maxdim, 1, nd);
            n = min(nx, prod(sz)); % Number of values in list
            idx = randperm(prod(sz));
            idx = idx(1:n)';
            
            X = tensor(@zeros,sz);
            X(idx) = 1;
            testCase.verifyEqual(X(idx, 'extract'), ones(n,1));
            testCase.verifyEqual(nnz(X), n);
            
            X = tensor(@zeros, sz);
            sub = tt_ind2sub(sz, idx);
            X(sub) = 1;
            testCase.verifyEqual(X(idx, 'extract'), ones(n,1));
            testCase.verifyEqual(nnz(X), n);
        end
        
        function SubsAsgnListArray(testCase, nd, maxdim, nx)
            sz = randi(maxdim, 1, nd);
            n = min(nx, prod(sz)); % Number of values in list
            idx = randperm(prod(sz));
            idx = idx(1:n)';
            
            X = tensor(@zeros,sz);
            X(idx) = (1:n)';
            testCase.verifyEqual(X(idx, 'extract'), (1:n)');
            testCase.verifyEqual(nnz(X), n);
            
            X = tensor(@zeros, sz);
            sub = tt_ind2sub(sz, idx);
            X(sub) = (1:n)';
            testCase.verifyEqual(X(idx, 'extract'), (1:n)');
            testCase.verifyEqual(nnz(X), n);
        end
        
        function SubsAsgnListArrayToEmpty(testCase, nd, maxdim, nx)
            sz = randi(maxdim, 1, nd);
            n = min(nx, prod(sz)); % Number of values in list
            idx = randperm(prod(sz));
            idx = idx(1:n)';
            
            % In the case that n = 1, we run into strange problems because
            % there is no way to differentiate between a 'list' and a
            % subscript when expanding an emtpy tensor. Maybe this can be
            % fixed by recognizing that it's a row rather than a column array?
            if (n > 1)
                X = tensor;
                sub = tt_ind2sub(sz, idx);
                X(sub) = (1:n)';
                newsz = max(sub,[],1);
                newidx = tt_sub2ind(newsz, sub);
                testCase.verifyEqual(X(newidx, 'extract'), (1:n)');
                testCase.verifyEqual(nnz(X), n);
            end
        end
        
%         function SubsAsgnListArrayError(testCase)
%             % Check for a linear index assignment that tries to expand the
%             % tensor. The problem is that it's not clear *how* to expand the
%             % tensor, so this should fail.
%             X = tensor(@zeros, [2 2 2]); %#ok<NASGU>
%             testCase.verifyError(@() eval('X(9) = 1;'), 'TTB:BadIndex');
%         end
        
        function SubsAsgnRange(testCase, nd, maxdim)
            sz = randi(maxdim, 1, nd);
            X = tensor(@zeros, sz);             %#ok<NASGU>
            iscolon = rand(nd,1) > 0.3; % Pick modes for range
            iscolon(randi(nd)) = 1; % Make sure there's at least one!
            newsz = sz;
            substrs = cell(nd,1);
            for k = 1:nd
                if iscolon(k)
                    if sz(k) == 1
                        substrs{k} = ':';
                    else
                        lower = randi(sz(k)-1);
                        upper = lower+randi(sz(k)-lower);
                        substrs{k} = sprintf('%d:%d', lower, upper);
                        newsz(k) = upper - lower + 1;
                    end
                else
                    substrs{k} = num2str(randi(sz(k)));
                end
            end
            substr = strjoin(substrs,',');
            newsz = newsz(iscolon);
            estr = ['X(' substr ') = tensor(@ones,newsz);'];
            eval(estr);
            testCase.verifyEqual(eval(['X(' substr ')']), tensor(@ones,newsz));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- FIND ---
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Find(testCase, szs)
            X = tenrand(szs);
            subs = find(X > 0.5);
            idx = find(X.data > 0.5);
            % This next step is a workaround since MATLAB's built-in function
            % returns a 0x1 array if X.data is nonempty but all zeros.
            if isempty(idx)
                idx = [];
            end
            testCase.verifyEqual(tt_sub2ind(size(X),subs), idx);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- TTM: Tensor Times Matrix ---
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function TtmOrder(testCase)
            % X x_i Mi x_j Mj = X x_j Mj x_i Mi
            nd = 4;
            sz = [4 3 2 4];
            X = tensor(@rand, sz);
            for i = 1:nd
                Mi = rand(5, sz(i));
                for j = (i+1):nd
                    Mj = rand(6, sz(j));
                    T1 = ttm( ttm(X,Mj,j), Mi, i);
                    T2 = ttm( ttm(X,Mi,i), Mj, j);
                    testCase.verifyEqual(T1.data, T2.data, 'RelTol', 1e-15);
                end
            end
        end
        
        function TtmSameMode(testCase)
            %X x1 A x1 B = X x1 BA
            nd = 4;
            sz = [4 3 2 4];
            X = tensor(@rand, sz);
            for i = 1:nd
                A = rand(5, sz(i));
                B = rand(6, 5);
                T1 = ttm( ttm(X,A,i), B, i);
                T2 = ttm( X, B*A, i);
                testCase.verifyEqual(T1.data, T2.data, 'RelTol', 1e-15);
            end
        end
        
        function Ttm(testCase)
            X = tensor(rand(5,3,4,2));
            A = rand(4,5); 
            B = rand(4,3); 
            C = rand(3,4); 
            D = rand(3,2);
            
            Y = ttm(X, A, 1); %<-- computes X times A in mode-1
            testCase.verifyEqual(size(Y), [4 3 4 2]);
            Z = A * reshape(X.data, 5, []);
            testCase.verifyEqual(Y.data(:), Z(:), 'RelTol', 1e-15);            
            Y2 = ttm(X, {A,B,C,D}, 1); %<-- same as above
            testCase.verifyEqual(Y.data, Y2.data, 'RelTol', 1e-15);            
            Y2 = ttm(X, A', 1, 't');   %<-- same as above
            testCase.verifyEqual(Y.data, Y2.data, 'RelTol', 1e-15);

            Y = ttm(X, {A,B,C,D}, [1 2 3 4]); %<-- 4-way multiply
            testCase.verifyEqual(size(Y), [4 4 3 3]);
            Y2 = ttm(X, {D,C,B,A}, [4 3 2 1]); %<-- same as above
            testCase.verifyEqual(Y.data, Y2.data, 'RelTol', 1e-15);
            Y2 = ttm(X, {A,B,C,D});
            testCase.verifyEqual(Y.data, Y2.data, 'RelTol', 1e-15);            
            Y2 = ttm(X, {A',B',C',D'}, 't');   %<-- same as above
            testCase.verifyEqual(Y.data, Y2.data, 'RelTol', 1e-15);            

            Y = ttm(X, {C,D}, [3 4]); %<-- X times C in mode-3 & D in mode-4
            testCase.verifyEqual(size(Y), [5 3 3 3]);         
            Y2 = ttm(X, {A,B,C,D}, [3 4]); %<-- same as above
            testCase.verifyEqual(Y.data, Y2.data, 'RelTol', 1e-15);            

            Y = ttm(X, {A,B,D}, [1 2 4]);   %<-- 3-way multiply
            testCase.verifyEqual(size(Y), [4 4 4 3]);         
            Y2 = ttm(X, {A,B,C,D}, [1 2 4]); %<-- same as above
            testCase.verifyEqual(Y.data, Y2.data, 'RelTol', 1e-15);            
            Y2 = ttm(X, {A,B,D}, -3);        %<-- same as above
            testCase.verifyEqual(Y.data, Y2.data, 'RelTol', 1e-15);
            Y2 = ttm(X, {A,B,C,D}, -3);      %<-- same as above
            testCase.verifyEqual(Y.data, Y2.data, 'RelTol', 1e-15);                         
        end
           
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- Symmetrize and Testing Symmetry ---
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function Symmetrize(testCase, nd, maxdim)
            n = randi(maxdim);
            sz = n * ones(1, nd);
            X = tensor(@rand, sz);
            S = symmetrize(X);
            if nd > 1 && prod(sz) > 1
                testCase.verifyFalse(issymmetric(X));
            end
            testCase.verifyTrue(issymmetric(S));
            sz = [4 3 2];
            X = tensor(@rand, sz);
            testCase.verifyError(@() symmetrize(X),'TTB:Tensor:BadModes');
        end
        
        function AllSymmetric(testCase,symver,issymver)
            m = 3;
            n = 4;
            sz = n * ones(1,m);
            X = tensor(rand(sz), sz);
            Y = symmetrize(X,1:m,symver);
            testCase.verifyTrue(issymmetric(Y,1:m,issymver));
        end
        
        function GroupedSymmetries(testCase,symver,issymver)
            X = tensor(rand(4,3,3,4));
            Y0 = symmetrize(X,{ [1 4], [2 3]},symver);
            testCase.verifyTrue(issymmetric(Y0,{[1 4], [2,3]},issymver));
        end
        
        function TestIsSymmetric(testCase,issymver)
            % Make sure it doesn't say everything is symmetric!
            m = 3;
            n = 4;
            sz = n * ones(1,m);
            X = tensor(rand(sz),sz);
            testCase.verifyFalse(issymmetric(X,1:m,issymver));
            
            sz = [4 3 2];
            X = tensor(rand(sz),sz);
            testCase.verifyFalse(issymmetric(X,1:m,issymver));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- TenFun ---
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        function Tenfun(testCase, nd, maxdim, gen)
            sz = randi(maxdim, 1, nd);
            X = tensor(gen, sz);
            Z = tenfun(@(x) x + 1, X);            
            testCase.verifyEqual(Z, X+1);
            
            Z = tenfun(@eq, X, 1);
            testCase.verifyEqual(Z, tensor(X.data == 1, sz));
            
            Y = tensor(gen, sz);
            Z = tenfun(@plus, X, Y);
            testCase.verifyEqual(Z, X+Y);
            
            W = tensor(gen,sz);
            Z = tenfun(@max, W, X, Y);
            T = max([X.data(:), Y.data(:), Z.data(:)], [], 2);
            testCase.verifyEqual(Z, tensor(T, sz));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- Collapse/SCALE ---
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function CollapseScale(testCase)
            X = tenrand([4 4 4]);
            Y = collapse(X, [2,3]);
            testCase.verifyEqual(Y.data, sum(reshape(X.data, [4 16]), 2));
            testCase.verifyEqual(size(Y), 4);
            Z = scale(X, 1./Y, 1);
            testCase.verifyEqual(double(collapse(Z, [2 3])), ones(4,1), 'RelTol', 1e-15);
            Y = collapse(X, [1 2], @max);
            testCase.verifyEqual(Y.data, max(reshape(X.data,[16 4]))');
            testCase.verifyEqual(size(Y), 4);
            Z = scale(X, 1./Y, 3);
            testCase.verifyEqual(double(collapse(Z, [1 2], @max)), ones(4,1), 'RelTol', 1e-15);
            X = tenones([3,4,5]);
            S = 10 * [1:5]'; 
            Y = scale(X,S,3);
            testCase.verifyEqual(double(scale(Y,1./S,3)), X.data, 'RelTol', 1e-15);
            S = tensor(10 * [1:5]',5); 
            Y = scale(X,S,3);
            testCase.verifyEqual(double(scale(Y,1./S,3)), X.data, 'RelTol', 1e-15);
            S = tensor(1:12,[3 4]); 
            Y = scale(X,S,[1 2]);
            testCase.verifyEqual(double(scale(Y,1./S,[1 2])), X.data, 'RelTol', 1e-15);
            S = tensor(1:12,[3 4]); 
            Y = scale(X,S,-3);
            testCase.verifyEqual(double(scale(Y,1./S,-3)), X.data, 'RelTol', 1e-15);
            S = tensor(1:60,[3 4 5]); 
            Y = scale(X,S,1:3);            
            testCase.verifyEqual(double(scale(Y,1./S,1:3)), X.data, 'RelTol', 1e-15);
        end
        
        function CollapseSum(testCase, szs)
            X = tenrand(szs);
            Y = collapse(X,2:ndims(X));
            Z = X.data;
            for j = 2:ndims(X)
                Z = sum(Z,j);
            end
            testCase.verifyEqual(double(Y),Z,'RelTol',1e-14);
        end
        
        function CollapseMax(testCase, szs)
            X = tenrand(szs);
            Y = collapse(X,1:ndims(X)-1,@max);
            Z = X.data;
            for j = 1:ndims(X)-1
                Z = max(Z,[],j);
            end
            if isempty(X.data)
                testCase.verifyEqual(double(Y),[],'RelTol',1e-14);
            else
                zsz = [size(X,ndims(X)) 1];
                testCase.verifyEqual(double(Y),reshape(Z, zsz),'RelTol',1e-14);
            end
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- Contract ---
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Contract(testCase)
            X = tensor(rand(4,3,2));
            Y = tensor(rand(3,2,4));
            Z1 = ttt(X,Y,1,3); %<-- Normal tensor multiplication
            Z2 = contract(ttt(X,Y),1,6); %<-- Outer product + contract
            testCase.verifyEqual(norm(Z1-Z2), 0, 'AbsTol', 1e-15);
        end
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- End ---
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function End(testCase)
            X = tensor(1:24, [4 3 2]);
            testCase.verifyEqual(X(end,1,1), X(4,1,1));
            testCase.verifyEqual(X(end,end,1), X(4,3,1));
            testCase.verifyEqual(X(end,end,end), X(4,3,2));
            testCase.verifyEqual(X(end), X(24));  
            testCase.verifyEqual(size(X(end,:,:)),[3 2]);
            testCase.verifyEqual(size(X(2:end,:,:)),[3 3 2]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- InnerProd ---
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Innerprod(testCase, nd, maxdim, gen)            
            sz = randi(maxdim, 1, nd);
            X = tensor(gen, sz);
            Y = tensor(gen, sz);
            Z = tensor(gen, sz+1);
            testCase.verifyEqual(innerprod(X, Y), dot(X.data(:), Y.data(:)), 'RelTol', 1e-14);
            testCase.verifyError(@() innerprod(X, Z), 'TTB:UnequalSize')                    
        end
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- MTTKRP ---
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Mttkrp(testCase)
            X = tensor(rand(2,3,4));
            A = rand(2,6);
            B = rand(3,6);
            C = rand(4,6);
            Y = mttkrp(X, {A,B,C}, 3);
            testCase.verifyEqual(size(Y), [4 6]);
            Y2 = mttkrp(X, {A,B,[]}, 3);
            testCase.verifyEqual(size(Y2), [4 6]);
            Z = double(tenmat(X,3))*khatrirao(B,A);
            testCase.verifyEqual(Y, Z, 'RelTol', 1e-14);
            testCase.verifyEqual(Y2, Z, 'RelTol', 1e-14);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- NVECS ---
        % This breaks in nd = 1 and sometimes if the dimension is 1.
        % Perhaps something to fix in a future version.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Nvecs(testCase, nd, maxdim)
            if nd > 1 && maxdim > 1
                sz = randi(maxdim, 1, nd);
                sz = max(sz,2);
                X = tensor(@rand, sz);
                n = randi(nd);
                r = randi(sz(n));
                U = nvecs(X, n, r);
                [U1,~,~] = svds(double(tenmat(X,n)),r);
                % Greedy sort
                for i=1:size(U1,2)
                    [~, j] = max(abs(U(:,i)'*U1(:,i:end)));
                    if (j ~= 1)
                        % Swap to remaining closest match.
                        U1(:, [i i+j-1]) = U1(:, [i+j-1 i]);
                    end
                    if (U(:,i)'*U1(:,i)<0)
                        % Fix direction.
                        U1(:,i) = -U1(:,i);
                    end
                end
                testCase.verifyEqual(U, U1, 'RelTol', 1e-8);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- PERMUTE ---
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Permute(testCase)
            X = tensor(rand(3,2,4));
            Y = permute(X, [1 3 2]);
            testCase.verifyEqual(size(Y), [3 4 2]);
            testCase.verifyEqual(Y.data, permute(X.data, [1 3 2]));
        end        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- RESHAPE ---
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Reshape(testCase)
            X = tensor(@rand, [4 3 2]);
            Y = reshape(X, [2 3 4]);
            testCase.verifyEqual(size(Y), [2 3 4]);
            testCase.verifyEqual(X.data(:), Y.data(:));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- SQUEEZE ---
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Squeeze(testCase)
            X = tensor(@rand, [2 1 3]);
            Y = squeeze(X);
            testCase.verifyEqual(size(Y), [2 3]);
            testCase.verifyClass(Y, 'tensor');
            X = tensor(@rand, [1 1]);
            Y = squeeze(X);
            testCase.verifyEqual(size(Y), [1 1]);
            testCase.verifyClass(Y, 'double');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- TTSV/TTV ---
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Ttv(testCase)
           X = tensor(@rand, [5,3,4,2]);
           A = rand(5,1); 
           B = rand(3,1); 
           C = rand(4,1); 
           D = rand(2,1);
           
           Y = ttv(X, A, 1); %<-- X times A in mode 1
           testCase.verifyEqual(size(Y), [3 4 2]);
           testCase.verifyClass(Y, 'tensor');
           Yalt = ttv(X, {A,B,C,D}, 1); %<-- same as above
           testCase.verifyEqual(Y.data, Yalt.data, 'RelTol', 1e-15);

           Y = ttv(X, {A,B,C,D}, [1 2 3 4]); %<-- All-mode multiply
           testCase.verifyEqual(size(Y), [1 1]);
           testCase.verifyClass(Y, 'double');
           Yalt = ttv(X, {D,C,B,A}, [4 3 2 1]); %<-- same as above
           testCase.verifyEqual(Y, Yalt);           
           Yalt = ttv(X, {A,B,C,D}); %<-- same as above
           testCase.verifyEqual(Y, Yalt);           
           
           Y = ttv(X, {C,D}, [3 4]); %<-- X times C in mode-3 & D in mode-4
           testCase.verifyEqual(size(Y), [5 3]);
           testCase.verifyClass(Y, 'tensor');
           Yalt = ttv(X, {A,B,C,D}, [3 4]); %<-- same as above
           testCase.verifyEqual(double(Y), double(Yalt), 'RelTol', 1e-15);

           Y = ttv(X, {A,B,D}, [1 2 4]); %<-- 3-way multiplication
           testCase.verifyEqual(size(Y), 4);
           testCase.verifyClass(Y, 'tensor');
           Yalt = ttv(X, {A,B,C,D}, [1 2 4]); %<-- same as above
           testCase.verifyEqual(double(Y), double(Yalt), 'RelTol', 1e-15);           
           Yalt = ttv(X, {A,B,D}, -3); %<-- same as above
           testCase.verifyEqual(double(Y), double(Yalt), 'RelTol', 1e-15);           
           Yalt = ttv(X, {A,B,C,D}, -3); %<-- same as above
           testCase.verifyEqual(double(Y), double(Yalt), 'RelTol', 1e-15);                      
        end
        
        
    end
    
    
    
end