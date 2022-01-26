% Testing importing tensors using import_data
classdef Test_ImportData < matlab.unittest.TestCase
    methods (Test)

        function Full(testCase)
            x = import_data('sptensor_small.tns');
            y = sptensor([1 1 1;2 2 2;3 3 3],[1 2 3]',[3 3 3]);
            testCase.verifyEqual(full(x), full(y));
        end
                    
    end
end