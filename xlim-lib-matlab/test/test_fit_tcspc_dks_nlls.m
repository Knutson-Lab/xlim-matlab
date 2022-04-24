% Test curves
irf = sim_tcspc_irfs(0.047,zeros(256,1),[20000, 0.8, 0.25, 0, 0],1);
trans = sim_tcspc_dks(0.047,irf,zeros(256,1),[20000, 2, 10000, 1, 0, 0, 0, 0],2,1);
wgts = trans;
wgts(trans == 0) = 1;
wgts(trans > 0) = 1 ./ trans(trans > 0);
%% Test 1: Proper lifetime extraction
results = fit_tcspc_dks_nlls(0.047,trans,irf,zeros(256,1),1,256,...
    [0.1, 3, 0.05, 1.5, 0, 0, 0, 0],vertcat(false(4,1),true(4,1)),...
    false(8,1),2,0,wgts,[],[],"EngineOptions",optimoptions("lsqnonlin",...
    "SpecifyObjectiveGradient",false));