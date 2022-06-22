classdef tcspcdata_test < matlab.unittest.TestCase
    methods(Test)
        function test_NumberOfBins(testCase)
            testdata = TcspcData("Data",ones(100,1),"BinWidth", 0.5, "PolarizationAngle", 44.7);
            testCase.assertEqual(testdata.NumberOfBins, 100);
        end
        function test_Data(testCase)
            testdata = TcspcData("Data",ones(100,1),"BinWidth", 0.5, "PolarizationAngle", 44.7);
            testCase.assertEqual(testdata.Data, ones(100,1));
        end
        function test_BinWidthError(testCase)
            testCase.verifyError(@()TcspcData("Data",ones(1,100),"BinWidth","a"),'MATLAB:validators:mustBeNumericOrLogical')
        end
    end
end