suite = testsuite("tcspcdata_test");

import matlab.unittest.TestRunner
runner = TestRunner.withTextOutput("OutputDetail",3);

import matlab.unittest.plugins.TestReportPlugin
pdfFile = "testreport.pdf";
p1 = TestReportPlugin.producingPDF(pdfFile);
runner.addPlugin(p1)

import matlab.unittest.plugins.XMLPlugin
xmlFile = "junittestresults.xml";
p2 = XMLPlugin.producingJUnitFormat(xmlFile);
runner.addPlugin(p2)

results = runner.run(suite);
