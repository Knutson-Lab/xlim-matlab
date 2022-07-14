classdef scalarparameterbounds_test < matlab.unittest.TestCase
    methods(Test)
        function test_boundValuesError(testCase)
            testCase.assertError(@()ScalarParameterBounds([1,5,3],[4,2,6]),'ScalarParameterBounds:constructor:ValueOfBounds')
        end  
        function test_boundSizeError(testCase)
            testCase.assertError(@()ScalarParameterBounds([1,2,3],[4,6]),'MATLAB:sizeDimensionsMustMatch')
        end 
        function test_multiGetLowerBound(testCase)
            bounds=ScalarParameterBounds([1,2,3],5);
            lb=bounds.multiGet("LowerBound");
            testCase.assertEqual(lb(2),2);
        end  
        function test_multiGetUpperBound(testCase)
            bounds=ScalarParameterBounds([],[5,6,7,8]);
            ub=bounds.multiGet("UpperBound");
            testCase.assertEqual(ub,[5;6;7;8]);
        end 
    end
end

