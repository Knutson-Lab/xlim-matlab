function [results,spa_results] = fit_tcspc_dks_nlls_subroutine...
    (tcal,ntcal,trans,ndata,ntrans,...
    irf,ntransirf,bck,ntransbck,fit_start,nfitstart,fit_end,...
    nfitend,param,nparam,ntau,ntheta,fxparam,gparam,spa_indx,...
    spa_low,spa_high,spa_size,lb,nlb,ub,nub,wgts,wgtflag,...
    psi,npsi,taur,ntaur,options,convflag,itpType)
    %FIT_TCSPC_DKS_NLLS_SUBROUTINE Drives the non-linear least squares
    %procedure

    %% Prepare engine arguments - Function handle
    % DATA - Structure that contains all relevant parameters for fit
    %%% TCAL - Bin width of data
    %%% NTCAL - Number of bin widths
    %%% TRANS - Data to compare to model
    %%% NDATA - Number of bins in data
    %%% NTRANS - Number of curves in model
    %%% IRF - IRF data to use in convolution
    %%% ITPIRF - Interpolated IRFs for color shifting
    itpirf = interp_irf(tcal,ntcal,irf,ndata,ntransirf,itpType);
    %%% NTRANSIRF - Number of IRFs
    %%% BCK - Background data to use for artifact control
    %%% NTRANSBCK - Number of backgrounds
    %%% WGT - Square root of weights to use when provided
    switch wgtflag
        case "Gaussian - Data"
            sqwgts = sqrt(gauss_wgts_data(trans));
        case "Custom"
            sqwgts = sqrt(wgts);
        otherwise
            sqwgts = wgts;
    end
    %%% NDATAWGT - Number of data points in weights
    %%% NTRANSWGT - Number of curves in weights
    [ndatawgt,ntranswgt] = size(sqwgts);
    %%% WGTFLAG - Flag for dynamic computation of weights
    %%% FIT_START - Where to start fitting
    %%% NFITSTART - Number of fit start values
    %%% FIT_END - Where to end fitting
    %%% NFITEND - Number of fit end values
    %%% PMAIN - Parameter vector for entire kernel
    %%% MPOINT - Pointer matrix of values for aligning PMAIN
    [pmain,mpoint,lbmain,ubmain] = get_pmain_mpoint_nlls...
        (param,nparam,fxparam,gparam,lb,nlb,ub,nub,ntrans);
    %%% FREEIDX - Free, unique indices in mpoint for indexing
    freeidx = unique(mpoint(mpoint > 0));
    %%% NPARAM - Number of parameters in model
    %%% NTAU - Number of lifetimes in model
    %%% NTHETA - Number of anisotropies in model
    %%% FXPARAM - Fixed parameters
    %%% GPARAM - Global parameters
    %%% TAUR - Reference lifetimes to use for standard based convolution
    %%% NTAUR - Number of reference lifetimes
    %%% PSI - Polarization angle at which data were acquired
    %%% NPSI - Number of polarization angles at which data were acquired
    %%% CONVFLAG - Flag indicating which convolution routine to use
    data = struct("tcal",tcal,"ntcal",ntcal,"trans",trans,"ndata",ndata,...
        "ntrans",ntrans,"irf",irf,"itpirf",itpirf,"ntransirf",ntransirf,...
        "bck",bck,"ntransbck",ntransbck,"wgt",sqwgts,"ndatawgt",ndatawgt,...
        "ntranswgt",ntranswgt,"wgtflag",wgtflag,"fit_start",fit_start,...
        "nfitstart",nfitstart,"fit_end",fit_end,"nfitend",nfitend,...
        "pmain",pmain,"mpoint",mpoint,"freeidx",freeidx,"nparam",nparam,...
        "ntau",ntau,"ntheta",ntheta,"fxparam",fxparam,"gparam",gparam,...
        "taur",taur,"ntaur",ntaur,"psi",psi,"npsi",npsi,"convflag",convflag);
    %%% F - Weighted residuals buffer for gated curves
    [F,J] = prep_engine_buffers(fit_start,nfitstart,fit_end,nfitend,...
        ntrans,numel(pmain));
    fun = @(x) fit_tcspc_dks_lsqnonlin_objfun(x,data,F,J);
    %% Prepare engine arguments - Other variables
    x0 = pmain(freeidx);
    if nlb > 0
        lb_engine = lbmain(freeidx);
    else
        lb_engine = lb;
    end

    if nub > 0
        ub_engine = ubmain(freeidx);
    else
        ub_engine = ub;
    end
    %% Run engine for main routine
    [x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin...
        (fun,x0,lb_engine,ub_engine,options);
    %% Organize outputs 
    results = fit_tcspc_dks_lsqnonlin_output(x,resnorm,residual,exitflag,...
        output,lambda,jacobian,data);
    spaflag = ~isempty(spa_indx);
    if spaflag && nargout > 1
        %% Generate SPA grid
        spagrids = cellfun(@(x,y,z) linspace(x,y,z),spa_low,spa_high,...
            spa_size,"UniformOutput",true);
        nspadims = numel(spa_low);
        [spa_coord{1:nspadims}] = ndgrid(spagrids{:});
        %% Initialize SPA results
        sparesultsgrid = cellfun(@(x) 1:x,spa_size,"UniformOutput",true);
        spa_results(sparesultsgrid{:}) = results;
        nevals = numel(spa_results);
        % Fix all parameter values in SPA range
        fxparam(spa_indx{:}) = true;
        for i = 1:nevals
            % Return parameters to fix
            spaparamvals = cellfun(@(x) x(i),spa_coord);
            % Assign to param & fix
            param(spa_indx{:}) = spaparamvals;
            % Recursively call this routine
            spa_results(i) = fit_tcspc_dks_nlls_subroutine...
                (tcal,ntcal,trans,ndata,ntrans,...
                irf,ntransirf,bck,ntransbck,fit_start,nfitstart,fit_end,...
                nfitend,param,nparam,ntau,ntheta,fxparam,gparam,spa_indx,...
                spa_low,spa_high,spa_size,lb,nlb,ub,nub,wgts,wgtflag,...
                psi,npsi,taur,ntaur,options,convflag,itpType);
        end
    end
end

function wgts = gauss_wgts_data(trans)
    % Subroutine that gets weights of 1 / N
    wgts = trans;
    nz = wgts > 0;
    wgts(nz) = 1 ./ wgts(nz);
    wgts(~nz) = 1;
end

function itpirf = interp_irf(tcal,ntcal,irf,ndata,ntransirf,itpType)
    % Interpolates IRF
    %% Allocate cell array to store cfit objects
    itpirf = cell(ntransirf,1);
    itpirf = fill_itpirf(itpirf,tcal,ntcal,irf,ndata,ntransirf,itpType);
end

function itpirf = fill_itpirf(itpirf,tcal,ntcal,irf,ndata,ntransirf,itpType)
    if ntcal == 1
        tax = (tcal / 2):tcal:(tcal * ndata);
    end
    for i = 1:ntransirf
        if ntcal > 1
            tcal_i = tcal(i);
            tax = (tcal_i / 2):tcal_i:(tcal_i * nbins);
        end
        itpirf{i} = fit(tax(:),irf(:,i),itpType);
    end
end

function [pmain,mpoint,lbmain,ubmain] = ...
    get_pmain_mpoint_nlls(param,nparam,fxparam,gparam,lb,nlb,ub,nub,ncurves)
    % Subroutine to compute parameter vector with pointer matrix
    %% Find mpoint first
    mpoint = get_mpoint(nparam,fxparam,gparam,ncurves);
    %% Return parameter vector
    Pm = repmat(param,1,ncurves);
    Mpa = abs(mpoint);
    mxind = max(Mpa,[],'all');
    pmain = zeros(mxind,1);
    pmain(Mpa) = Pm;
    if nlb > 0
        Pm = repmat(lb,1,ncurves);
        lbmain = zeros(mxind,1);
        lbmain(Mpa) = Pm;
    else
        lbmain = [];
    end
    if nub > 0
        Pm = repmat(ub,1,ncurves);
        ubmain = zeros(mxind,1);
        ubmain(Mpa) = Pm;
    else
        ubmain = [];
    end
    
end

function [B_resid,B_jac] = prep_engine_buffers...
    (fit_start,nfitstart,fit_end,nfitend,ncurves,ntotalparam)
    chanwindow = fit_end - fit_start + 1;
    if nfitstart == 1 && nfitend == 1
        ntotalbins = ncurves * chanwindow;
    else
        ntotalbins = sum(chanwindow);
    end
    B_resid = zeros(ntotalbins,1);
    B_jac = zeros(ntotalbins,ntotalparam);
end