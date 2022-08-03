classdef nonlinearleastsquaresparameter_test < matlab.unittest.TestCase
    methods(Test)
        function test_ValuesOutsideBoundsError(testCase)
            testscalarbounds = ScalarParameterBounds(1,10);
            testconstraints = BoundaryConstraints("ScalarBounds",testscalarbounds);
            testfitparam = FitParameter("Value",12);

            testCase.assertError(@()NonLinearLeastSquaresParameter(testfitparam,testconstraints),'OptimizationParameter:constructor:ValueOutsideBounds')
        end  
        function test_multiSetConstraintsError(testCase)

            testscalarbounds = ScalarParameterBounds(1,10);
            testconstraints = BoundaryConstraints("ScalarBounds",testscalarbounds);
            testfitparam = FitParameter("Value",5,"Fixed",false,"Linked",true);

            nlls2(2) = NonLinearLeastSquaresParameter;

            testCase.assertError(@()nlls2.multiSet("Value",[testfitparam,testfitparam],"Constraint",[testconstraints,testconstraints,testconstraints]),'NonLinearLeastSquaresParameter:multiSet:ConstarintsSize')
        end
        function test_multiGetConstraintsError(testCase)

            testscalarbounds = ScalarParameterBounds(1,10);
            testconstraints = BoundaryConstraints("ScalarBounds",testscalarbounds);
            testfitparam = FitParameter("Value",5,"Fixed",false,"Linked",true);

            nlls2(2) = NonLinearLeastSquaresParameter;
            nlls2=nlls2.multiSet("Value",[testfitparam,testfitparam],"Constraint",[testconstraints,testconstraints]);
            [vals,consts] = nlls2.multiGet("Value","Constraints");

            testCase.assertEqual(vals(1).Value,5);
        end
    end
end
