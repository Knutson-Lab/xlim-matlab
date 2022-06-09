function tests = tcspcdata_test
    tests = functiontests(localfunctions);
end

function test_NumberOfBins(testCase)
    testdata = TcspcData("Data",ones(100,1),"BinWidth", 0.5, "PolarizationAngle", 44.7);
    testCase.assertEqual(testdata.NumberOfBins, 100);
end

function test_Data(testCase)
    testdata = TcspcData("Data",ones(100,1),"BinWidth", 0.5, "PolarizationAngle", 44.7);
    testCase.assertEqual(testdata.Data, ones(100,1));
end
