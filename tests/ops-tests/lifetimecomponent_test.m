classdef lifetimecomponent_test < matlab.unittest.TestCase
    methods(Test)
        function test_WeightedLifetimeComponentEqual(testCase)
            testscalarbounds = ScalarParameterBounds(0,10);
            testconstraints = BoundaryConstraints("ScalarBounds",testscalarbounds);
            testfitparam1 = FitParameter("Value",1);
            testfitparam2 = FitParameter("Value",4,"Fixed",true,"Linked",false);

            testtau = NonLinearLeastSquaresParameter(testfitparam2,testconstraints);
            testamp = NonLinearLeastSquaresParameter(testfitparam1,testconstraints);
            testwlifetime = WeightedLifetimeComponent("Amplitude",testamp,"Lifetime",testtau);

            testCase.assertEqual(testwlifetime.Lifetime.Value.Value,4);
        end
    end
end

