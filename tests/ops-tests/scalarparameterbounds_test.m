classdef scalarparameterbounds_test < matlab.unittest.TestCase
    methods(Test)
        function test_FixedBoolean(testCase)
            bounds=ScalarParameterBounds([1,2,3],[4,5,6]);
            testCase.assertEqual(bounds.LowerBound(2),2);
        end    
    end
end