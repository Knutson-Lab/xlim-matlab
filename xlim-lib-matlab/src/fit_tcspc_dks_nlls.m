function [results,spa_results] = fit_tcspc_dks_nlls...
    (tcal,trans,irf,bck,fit_start,fit_end,param,fxparam,gparam,ntau,...
    ntheta,spa_indx,spa_low,spa_high,spa_size,args)
    %FIT_TCSPC_DKS_NLLS Fits curves according to a non-linear least squares
    %multidimensional minimization routine
    
    arguments
        tcal (:,1) {mustBePositive,mustBeFinite}
        trans (:,:) {mustBeNonnegative,mustBeInteger}
        irf (:,:) {mustBeNonnegative,mustBeInteger}
        bck (:,:) {mustBeNonnegative,mustBeInteger}
        fit_start (:,1) {mustBePositive,mustBeInteger}
        fit_end (:,1) {mustBePositive,mustBeInteger}
        param (:,1) {mustBeNumeric,mustBeFinite}
        fxparam (:,1) logical
        gparam (:,1) logical
        ntau (1,1) {mustBeNonnegative,mustBeInteger}
        ntheta (1,1) {mustBeNonnegative,mustBeInteger}
    end

    arguments (Repeating)
        spa_indx (1,1) {mustBePositive,mustBeInteger}
        spa_low (1,1) {mustBeNumeric,mustBeFinite}
        spa_high (1,1) {mustBeNumeric,mustBeFinite}
        spa_size (1,1) {mustBePositive,mustBeInteger}
    end

    arguments
        args.LowerBounds (:,1) {mustBeNumeric} = []
        args.UpperBounds (:,1) {mustBeNumeric} = []
        args.Weights (:,:) {mustBeNonnegative,mustBeFinite}
        args.WeightMethod {mustBeMember(args.WeightMethod,["Gaussian - Data",...
            "Custom"])} = "Gaussian - Data"
        args.PolarizationAngle (:,1) {mustBeInRange...
            (args.PolarizationAngle,0,90)} = 44.7
        args.ReferenceLifetime (:,1) {mustBePositive,mustBeFinite}
        args.EngineOptions (1,1) optim.options.Lsqnonlin = optimoptions("lsqnonlin")
        args.ConvolutionMethod {mustBeMember(args.ConvolutionMethod,...
            "FFT")} = "FFT"
        args.InterpolationOption {mustBeMember(args.InterpolationOption,...
            ["smoothingspline","cubicinterp"])} = "smoothingspline"
    end
    %% Validate required argument sizes - Curves
    ntcal = numel(tcal);
    [ndata,ntrans] = size(trans);
    [ndatairf,ntransirf] = size(irf);
    [ndatabck,ntransbck] = size(bck);
    nfitstart = numel(fit_start);
    nfitend = numel(fit_end);
    validate_req_arg_sz_curves(ntcal,ndata,ntrans,ndatairf,ntransirf,...
        ndatabck,ntransbck,nfitstart,nfitend);
    %% Validate argument bounds - Gates
    % FIT_START/FIT_END - fitstart < fitend <= ndata
    assert(all(fit_start < fit_end) && all(fit_end <= ndata))
    %% Validate required argument sizes - Parameters
    nparam = numel(param);
    % FXPARAM/GPARAM - nparam x 1
    nfxparam = numel(fxparam);
    ngparam = numel(gparam);
    assert(nfxparam == nparam)
    assert(ngparam == nparam);
    % NTAU/NTHETA - (2 * ntau) + (2 * ntheta) + 4 == nparam
    nexpparam = (2 * ntau) + (2 * ntheta) + 4;
    assert(nexpparam == nparam);
    %% Valdiate support plane block
    % SPA_INDX - Between 1 - nparam
    assert(all(cellfun(@(x) x <= nparam,spa_indx)));
    % SPA_LOW/SPA_HIGH - spa_low < spa_high
    assert(all(cellfun(@(x,y) x < y,spa_low,spa_high)))
    %% Validate optional arguments block - Bounds
    % Lower/upper bounds - Must be nparam x 1
    nlb = numel(args.LowerBounds);
    nub = numel(args.UpperBounds);
    if nlb > 0
        assert(nlb == nparam);
        assert(args.LowerBounds <= param);
    end
    if nub > 0
        assert(nub == nparam);
        assert(param <= args.UpperBounds);
    end
    %% Validate optional arguments block - Weights
    % wgts - When custom, ndata x 1, ndata x ntrans, 1 x 1, or 1 x ntrans
    if isfield(args,"Weights")
        [ndatawgts,ntranswgts] = size(args.Weights);
        assert(args.WeightMethod == "Custom");
        assert((ndatawgts == ndata || ndatawgts == 1) && ...
            (ntranswgts == ntrans || ntranswgts == 1))
    else
        assert(args.WeightMethod ~= "Custom");
        args.Weights = [];
    end
    %% Validate optional arguments block - Polarization
    % psi - ndata x 1, or 1 x 1
    npsi = numel(args.PolarizationAngle);
    assert(npsi == ntrans || npsi == 1);
    %% Validate optional arguments block - Reference lifetime
    % taur - ndata x 1, or 1 x 1
    if ~isfield(args,"ReferenceLifetime")
        args.ReferenceLifetime = [];
    end
    ntaur = numel(args.ReferenceLifetime);
    if ntaur > 0
        assert(ntaur == ntrans || ntaur == 1);
    end
    %% Validate optional arguments block - Engine
    % Cannot run if color shift value is free and engine has an objective
    % gradient
    fxq = fxparam(nparam - 3);
    assert(~(~fxq && args.EngineOptions.SpecifyObjectiveGradient))
    %% Pass to subroutine
    if nargout == 1
        [results] = ...
            fit_tcspc_dks_nlls_subroutine(tcal,ntcal,trans,ndata,ntrans,...
            irf,ntransirf,bck,ntransbck,fit_start,nfitstart,fit_end,...
            nfitend,param,nparam,ntau,ntheta,fxparam,gparam,spa_indx,...
            spa_low,spa_high,spa_size,args.LowerBounds,nlb,args.UpperBounds,...
            nub,args.Weights,args.WeightMethod,args.PolarizationAngle,npsi,...
            args.ReferenceLifetime,ntaur,args.EngineOptions,...
            args.ConvolutionMethod,args.InterpolationOption);
    else
        [results,spa_results] = ...
            fit_tcspc_dks_nlls_subroutine(tcal,ntcal,trans,ndata,ntrans,...
            irf,ntransirf,bck,ntransbck,fit_start,nfitstart,fit_end,...
            nfitend,param,nparam,ntau,ntheta,fxparam,gparam,spa_indx,...
            spa_low,spa_high,spa_size,args.LowerBounds,nlb,args.UpperBounds,...
            nub,args.Weights,args.WeightMethod,args.PolarizationAngle,npsi,...
            args.ReferenceLifetime,ntaur,args.EngineOptions,...
            args.ConvolutionMethod,args.InterpolationOption);
    end
    
end

function validate_req_arg_sz_curves(ntcal,ndata,ntrans,ndatairf,ntransirf,...
    ndatabck,ntransbck,nfitstart,nfitend)
    % TCAL - ntrans x 1, or 1 x 1
    assert(ntcal == ntrans || ntcal == 1)
    % IRF - ndata x 1, or ndata x ntrans
    assert(ndata == ndatairf && (ntransirf == ntrans || ntransirf == 1))
    % BCK - Same as IRF
    assert(ndata == ndatabck && (ntransbck == ntrans || ntransbck == 1))
    % FIT_START/FIT_END - ntrans x 1, or 1 x 1
    assert(nfitstart == ntrans || nfitstart == 1)
    assert(nfitend == ntrans || nfitend == 1)
end
