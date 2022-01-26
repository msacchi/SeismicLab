% Testing importing tensors using import_data
classdef Test_KtensorScore < matlab.unittest.TestCase
    methods (Test)

        function ScoreDifferentKtensors(testCase)
            K1 = ktensor([1; 1], [1, 2; 3, 4], [5, 6; 7, 8]);
            K2 = ktensor([1; 1], [0, 2; 3, 4], [5, 6; 7, 8]);
            testCase.verifyEqual(score(K1,K2), 0.95, 'RelTol', 1e-15);
        end
                    
        function ScoreSameKtensors(testCase)
            K1 = ktensor([1; 1], [1, 2; 3, 4], [5, 6; 7, 8]);
            testCase.verifyEqual(score(K1,K1), 1.0, 'RelTol', 1e-15);
        end
        
    end
end
