% Testing importing tensors using import_data
classdef Test_ExportData < matlab.unittest.TestCase
    methods (Test)

        function Full(testCase)
            x = import_data('ktensor_small.ktns');
            export_data(x,'ktensor_small_2.ktns');
            x2 = import_data('ktensor_small_2.ktns');
            testCase.verifyEqual(x, x2);
            if exist('ktensor_small_2.ktns', 'file')==2
                delete('ktensor_small_2.ktns');
            end
        end % Full
        
        function Formatting_Integer_Sptensor(testCase)
            X = sptensor([1 1 1;2 2 2; 4 3 2],[1;2;3],[4 3 2]);
            export_data(X,'X_default.sptensor');
            export_data(X,'X_int_data.sptensor','fmt_data','%d');
            X1 = import_data('X_default.sptensor');
            X2 = import_data('X_int_data.sptensor');
            testCase.verifyEqual(X1,X);
            testCase.verifyEqual(X2,X);
            testCase.verifyEqual(X1,X2);
            if exist('X_default.sptensor', 'file')==2
                delete('X_default.sptensor');
            end
            if exist('X_int_data.sptensor', 'file')==2
                delete('X_int_data.sptensor');
            end
        end % Formatting_Integer_Sptensor

            
        function Formatting_Integer_Ktensor(testCase)
            K = ktensor([1; 2], 4*ones(4,2), 5*ones(5,2), 3*ones(3,2));
            export_data(K,'K_default.ktensor'); 
            export_data(K,'K_int_data.ktensor','fmt_data','%d');
            export_data(K,'K_int_data_int_lambda.ktensor','fmt_data','%d','fmt_lambda','%d');
            K1 = import_data('K_default.ktensor');
            K2 = import_data('K_int_data.ktensor');
            K3 = import_data('K_int_data_int_lambda.ktensor');
            testCase.verifyEqual(K1,K);
            testCase.verifyEqual(K2,K);
            testCase.verifyEqual(K3,K);
            if exist('K_default.ktensor', 'file')==2
                delete('K_default.ktensor');
            end
            if exist('K_int_data.ktensor', 'file')==2
                delete('K_int_data.ktensor');
            end
            if exist('K_int_data_int_lambda.ktensor', 'file')==2
                delete('K_int_data_int_lambda.ktensor');
            end
        end % Formatting_Integer_Ktensor
    
        
    end % methods
end % classdef Test_ExportData
