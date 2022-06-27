classdef tcspcinstrument_test < matlab.unittest.TestCase
    methods(Test)
        function test_RefLifetimeError(testCase)
            testdata = TcspcData("Data",ones(100,1),"BinWidth", 0.5, "PolarizationAngle", 44.7);
            testCase.verifyError(@()TcspcInstrument("Data",testdata,"ReferenceLifetime",-2),'MATLAB:validators:mustBePositive');
        end
        function test_IrfData(testCase)
            testdata = TcspcData("Data",ones(100,1),"BinWidth", 0.5, "PolarizationAngle", 44.7);
            testirf = TcspcInstrument("Data",testdata);
            testCase.assertEqual(testirf.Data, testdata);
        end
    end
end