classdef tcspcinstrument_test < matlab.unittest.TestCase
    methods(Test)
        function test_RefLifetimeError(testCase)
            testdata = TcspcData("Data",ones(100,1),"BinWidth", 0.5, "PolarizationAngle", 44.7);
            testCase.verifyError(@()TcspcInstrument("Data",testdata,"ReferenceLifetime",-2),'MATLAB:validators:mustBePositive');
        end
        function test_IrfData(testCase)
            testdata = TcspcData("Data",ones(100,1),"BinWidth", 0.5, "PolarizationAngle", 44.7);
            testirf = TcspcInstrument("Data",testdata,"ReferenceLifetime",2);
            testCase.assertEqual(testirf.Data, testdata);
        end

        function test_multiSetDataError(testCase)
            data100(100) = TcspcData;
            data100 = data100.multiSet("Data",ones(100),"BinWidth",0.5*ones(100,1),"PolarizationAngle", 44.7);
            irf1000(1000) = TcspcInstrument;

            testCase.assertError(@()irf1000.multiSet("Data",data100,"ReferenceLifetime",2),'MATLAB:assertion:failed')
        end
        function test_multiGet_getIrfData(testCase)
            data1000(1000) = TcspcData;
            data1000=data1000.multiSet("Data",ones(1000),"BinWidth",0.5*ones(1000,1),"PolarizationAngle",44.7);
            irf1000(1000) = TcspcInstrument;
            irf1000=irf1000.multiSet("Data",data1000,"ReferenceLifetime",2);

            irfset=irf1000.multiGet("Data");
            dataset=irfset.multiGet("Data");
            testCase.assertEqual(dataset, ones(1000));
        end
    end
end