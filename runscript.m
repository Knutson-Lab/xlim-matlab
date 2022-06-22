suite = testsuite("xlim-ops-matlab/test");
addpath(genpath('xlim-ops-matlab'))

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

import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoberturaFormat
%sourceCodeFile = "xlim-ops-matlab/src/@TcspcData/TcspcData.m";
sourceCodeFolder = "xlim-ops-matlab/src";
reportFile = "cobertura.xml";
reportFormat = CoberturaFormat(reportFile);
p3 = CodeCoveragePlugin.forFolder(sourceCodeFolder,"Producing",reportFormat);
runner.addPlugin(p3)

results = runner.run(suite);
