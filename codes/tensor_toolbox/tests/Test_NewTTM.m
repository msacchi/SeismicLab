% Testing different versions of tensor/ttm for correctness and efficiency
classdef Test_NewTTM < matlab.unittest.TestCase
    
    
    properties (TestParameter)
        combo = struct('small3d', [3 10 50], 'small4d', [4 10 25], 'small5d', [5 5 10], 'large3d',[3 100 250]);
        ver = struct('old', 0, 'new', 1);
    end
    
    methods (Test)
        function Compare(testCase, combo)
            nd = combo(1);
            lsz = combo(2);
            usz = combo(3);
            rsz = usz - lsz;
            sz = lsz * ones(1, nd) + randi(rsz, 1, nd);
            X = tensor(@randn, sz);
            U = cell(nd,1);
            for n = 1:nd
                U{n} = randn(lsz + randi(rsz), sz(n));
            end
            for n = 1:nd
                Y1 = ttm(X,U{n},n,[],0);
                Y2 = ttm(X,U{n},n,[],1);
                testCase.verifyEqual(size(Y1), size(Y2));
                testCase.verifyEqual(Y1.data, Y2.data, 'AbsTol', 1e-12);
            end      
        end
        
        function Time(testCase, combo, ver)
            nd = combo(1);
            lsz = combo(2);
            usz = combo(3);
            sz = usz * ones(1, nd);
            X = tensor(@randn, sz);
            U = cell(nd,1);
            for n = 1:nd
                U{n} = randn(lsz, sz(n));
            end
            for n = 1:nd
                newsz = sz;
                newsz(n) = lsz;
                Y1 = ttm(X,U{n},n,[],ver);
                testCase.verifyEqual(size(Y1), newsz);
            end      
        end
        
    end
end
