classdef boundaryconstraints_test < matlab.unittest.TestCase

    methods(Test)
        function test_constructor(testCase)
            scalarbounds = ScalarParameterBounds(1,5);
            testbounds = BoundaryConstraints("ScalarBounds",scalarbounds,"LinearInequalityVector",zeros(2,1),...
                "LinearInequalityScalar",0,"LinearEqualityVector",ones(2,1),"LinearEqualityScalar",1);
            testCase.assertEqual(testbounds.ScalarBounds.LowerBound,1);
        end  
    end
end

