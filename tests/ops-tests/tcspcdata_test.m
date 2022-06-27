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
        function test_multiSetDataError(testCase)
            %correct format 
            %data1000(1000) = TcspcData
            %data1000=data1000.multiSet("Data",ones(1000),"BinWidth",0.5*ones(1000,1),"PolarizationAngle",44.7)
            data1000(1000) = TcspcData;
            testCase.assertError(@()data1000.multiSet("Data",ones(100),"BinWidth",0.5*ones(1000,1),"PolarizationAngle",44.7),'MATLAB:assertion:failed')
        end
        function test_multiSetBinWidthError(testCase)
            data1000(1000) = TcspcData;
            testCase.assertError(@()data1000.multiSet("Data",ones(1000),"BinWidth",0.5*ones(100,1),"PolarizationAngle",44.7*ones(1000,1)),'MATLAB:assertion:failed')
        end
        function test_multiSetPolAngError(testCase)
            data1000(1000) = TcspcData;
            testCase.assertError(@()data1000.multiSet("Data",ones(1000),"BinWidth",0.5*ones(1000,1),"PolarizationAngle",[44.7,44.7]),'MATLAB:assertion:failed')
        end
        function test_multiGet_getData(testCase)
            data1000(1000) = TcspcData;
            data1000=data1000.multiSet("Data",ones(1000),"BinWidth",0.5*ones(1000,1),"PolarizationAngle",44.7);
            dataset=data1000.multiGet("Data");
            testCase.assertEqual(dataset, ones(1000));
        end
    end
end