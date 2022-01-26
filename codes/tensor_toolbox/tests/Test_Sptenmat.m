% Testing importing tensors using import_data
classdef Test_Sptenmat < matlab.unittest.TestCase
    methods (Test)

        function SptenmatSptensorInput(testCase)
            % sparse tensor
            T = sptensor([1 1 1; 2 2 2; 3 3 3],[1; 2; 3],[3 3 3]);
            % sparse tensor as matrix
            M = sptenmat(T,1);
            % equivalent sparse matrix
            Msp = sparse([1 2 3],[1 5 9],[1 2 3], 3, 9);
            % test equality by converting both matrices to dense, double
            testCase.verifyEqual(full(double(M)), full(Msp));
        end
        
        function SptenmatInvalidInput(testCase)
            % dense tensor
            T = tenrand([2 3 4]);
            % cannot use dense tensor as input to sptenmat
            testCase.verifyError(@()sptenmat(T,1), 'sptenmat:InvalidInput');
        end
    end
end
