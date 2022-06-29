classdef fitparameter_test < matlab.unittest.TestCase
    methods(Test)
        function test_ValueError(testCase)
            testCase.verifyError(@()FitParameter("Value","a","Fixed",0,"Linked",0),'MATLAB:validators:mustBeNumeric')
        end
        function test_multiSetLinkedError(testCase)
            param100(100)=FitParameter;
            testCase.verifyError(@()param100.multiSet("Value",ones(100,1),"Fixed",1,"Linked",[0,1]),'MATLAB:assertion:failed')
        end
        function test_multiGet_getValue(testCase)
            param100(100)=FitParameter;
            param100=param100.multiSet("Value",ones(100,1),"Fixed",ones(1,100),"Linked",0);
            valueset=param100.multiGet("Value");
            testCase.assertEqual(valueset, ones(100,1));
        end
    end
end