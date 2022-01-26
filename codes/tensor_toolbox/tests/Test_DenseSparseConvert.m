% Testing conversions between tensor and sptensor
classdef Test_DenseSparseConvert < matlab.unittest.TestCase
    methods (Test)

        function Empty(testCase)
            x = sptensor;
            y = tensor;
            testCase.verifyEqual(x, sptensor(y));
            testCase.verifyEqual(y, tensor(x));
        end
            
        function Zero(testCase)
            x = sptensor([5 4 3]);
            y = tensor(@zeros, [5 4 3]);
            testCase.verifyEqual(x, sptensor(y));
            testCase.verifyEqual(y, tensor(x));
        end
        
        function ThreeWay(testCase)
            x = sptenrand([4 3 2], 0.4);
            y = tensor(x);
            testCase.verifyEqual(x, sptensor(y));
            
            y = tenrand([4 3 2]);
            x = sptensor(y);
            testCase.verifyEqual(y, tensor(x));            
        end
        
        function OneWay(testCase)
            x = sptenrand(10,0.4);
            y = tensor(x);
            testCase.verifyEqual(x, sptensor(y));
            
            y = tenrand([10]);
            x = sptensor(y);
            testCase.verifyEqual(y, tensor(x));            
        end
        
    end
end