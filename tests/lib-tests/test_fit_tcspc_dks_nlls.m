% % Test curves
% irf = sim_tcspc_irfs(0.047,zeros(256,1),[20000, 0.8, 0.25, 0, 0],1);
% trans = horzcat(sim_tcspc_dks(0.047,irf,zeros(256,1),[200000, 2, 100000, 1, 0, 0, 0, 0],2,10),sim_tcspc_dks(0.047,irf,zeros(256,1),[100000, 2, 200000, 1, 0, 0, 0, 0],2,10));
%% Test 1: Proper lifetime extraction
% results = fit_tcspc_dks_nlls(0.047,trans,irf,zeros(256,1),1,256,...
%     [0.1, 3, 0.05, 1.5, 0, 0, 0, 0],vertcat(false(4,1),true(4,1)),...
%     [false true false true false(1,4)],2,0,"EngineOptions",optimoptions("lsqnonlin",...
%     "SpecifyObjectiveGradient",true));
% 
% testCase.verifyError(@()fit_tcspc_dks_nlls(0.047,trans,irf,zeros(256,1),1,256,...
%     [0.1, 3, 0.05, 1.5, 0, 0, 0, 0],vertcat(false(4,1),true(4,1)),...
%     [false true false true false(1,4)],2,0,"EngineOptions",optimoptions("lsqnonlin",...
%     "SpecifyObjectiveGradient",true)),'MATLAB:validators:mustBeNumericOrLogical')

classdef test_fit_tcspc_dks_nlls < matlab.unittest.TestCase
    methods(Test)
        function test_firsttry(testCase)

            irf = sim_tcspc_irfs(0.047,zeros(256,1),[20000, 0.8, 0.25, 0, 0],1);
            trans = horzcat(sim_tcspc_dks(0.047,irf,zeros(256,1),[200000, 2, 100000, 1, 0, 0, 0, 0],2,10),...
            sim_tcspc_dks(0.047,irf,zeros(256,1),[100000, 2, 200000, 1, 0, 0, 0, 0],2,10));

            testCase.verifyError(@()fit_tcspc_dks_nlls(0.047,trans,irf,zeros(256,1),1,256,...
                [0.1, 3, 0.05, 1.5, 0, 0, 0, 0],vertcat(false(4,1),true(4,1)),...
                [false true false true false(1,4)],2,0,"EngineOptions",optimoptions("lsqnonlin",...
                "SpecifyObjectiveGradient",true)),'MATLAB:assertion:failed')            
        end
    end
end    

%% Test 2: SPA analysis
% [~,spa_results] = fit_tcspc_dks_nlls(0.047,trans,irf,zeros(256,1),1,256,...
%     [0.1, 3, 0.05, 1.5, 0, 0, 0, 0],vertcat(false(4,1),true(4,1)),...
%     [false false false false false(1,4)],2,0,2,1,3,20,...
%     "EngineOptions",optimoptions("lsqnonlin",...
%     "SpecifyObjectiveGradient",true));