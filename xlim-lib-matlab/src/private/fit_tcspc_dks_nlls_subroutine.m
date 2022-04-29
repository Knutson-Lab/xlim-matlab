function [results,spa_results] = fit_tcspc_dks_nlls_subroutine...
    (tcal,ntcal,trans,ndata,ntrans,...
    irf,ntransirf,bck,ntransbck,fit_start,nfitstart,fit_end,...
    nfitend,param,nparam,ntau,ntheta,fxparam,gparam,spa_indx,...
    spa_low,spa_high,spa_size,lb,nlb,ub,nub,wgts,wgtflag,...
    psi,npsi,taur,ntaur,options,convflag,itpType,parflag,typeflag)
    %FIT_TCSPC_DKS_NLLS_SUBROUTINE Drives the non-linear least squares
    %procedure

    %% Check workflow
    
    if typeflag == "Single"
        % Single curve-by-curve workflow
        [~,ntranswgt] = size(wgts);
        
        % Initialize output results structures for curves
        results(1:ntrans,1) = struct("curves",...
            struct("fitted",zeros(ndata,1),...
            "residuals",zeros(ndata,1),"param",zeros(nparam,1),...
            "sterr",zeros(nparam,1),"X2",0,"df",0),...
            "X2",0,"df",0,"residuals",zeros(ndata,1),"exitflag",0,...
            "output",struct([]),"lambda",struct([]),"jacobian",...
            zeros(ndata,numel(find(~fxparam))));
        if nargout > 1
            %% Initialize SPA results
            sparesultsgrid = cellfun(@(x) 1:x,spa_size,"UniformOutput",false);
            spa_results_l(sparesultsgrid{:}) = results(1);
            spa_results(1:ntrans,1) = struct("spagrid",spa_results_l);
        end

        if parflag
            nout = narargout;
            parfor i = 1:ntrans
                % Pass recursively with simultaneous flag
                tcal_l = tcal;
                if ntcal > 1
                    tcal_i = tcal_l(i);
                else
                    tcal_i = tcal_l;
                end
                ntcal_i = 1;
                trans_i = trans(:,i);
                irf_l = irf;
                if ntransirf > 1
                    irf_i = irf_l(:,i);
                else
                    irf_i = irf_l;
                end
                ntransirf_i = 1;
                bck_l = bck;
                if ntransbck > 1
                    bck_i = bck_l(:,i);
                else
                    bck_i = bck_l;
                end
                ntransbck_i = 1;
                fit_start_l = fit_start;
                if nfitstart > 1
                    fit_start_i = fit_start_l(i);
                else
                    fit_start_i = fit_start_l;
                end
                nfitstart_i = 1;
                fit_end_l = fit_end;
                if nfitend > 1
                    fit_end_i = fit_end_l(i);
                else
                    fit_end_i = fit_end_l;
                end
                nfitend_i = 1;
                wgts_l = wgts;
                if ntranswgt > 1
                    wgts_i = wgts_l(:,i);
                else
                    wgts_i = wgts_l;
                end
                psi_l = psi;
                if npsi > 1
                    psi_i = psi_l(i);
                else
                    psi_i = psi_l;
                end
                taur_l = taur;
                npsi_i = 1;
                if ntaur > 1
                    taur_i = taur_l(i);
                else
                    taur_i = taur_l;
                end
                ntaur_i = 1;
                if nout > 1
                    [results(i),spa_results(i)] = ...
                        fit_tcspc_dks_nlls_subroutine...
                        (tcal_i,ntcal_i,trans_i,ndata,1,...
                        irf_i,ntransirf_i,bck_i,ntransbck_i,fit_start_i,...
                        nfitstart_i,fit_end_i,nfitend_i,param,nparam,...
                        ntau,ntheta,fxparam,gparam,spa_indx,...
                        spa_low,spa_high,spa_size,lb,nlb,ub,nub,wgts_i,wgtflag,...
                        psi_i,npsi_i,taur_i,ntaur_i,options,convflag,...
                        itpType,false,"Simultaneous");  
                else
                    [results(i)] = ...
                        fit_tcspc_dks_nlls_subroutine...
                        (tcal_i,ntcal_i,trans_i,ndata,1,...
                        irf_i,ntransirf_i,bck_i,ntransbck_i,fit_start_i,...
                        nfitstart_i,fit_end_i,nfitend_i,param,nparam,...
                        ntau,ntheta,fxparam,gparam,spa_indx,...
                        spa_low,spa_high,spa_size,lb,nlb,ub,nub,wgts_i,wgtflag,...
                        psi_i,npsi_i,taur_i,ntaur_i,options,convflag,...
                        itpType,false,"Simultaneous");
                end
                
            end
        else
            for i = 1:ntrans
                % Pass recursively with simultaneous flag
                if ntcal > 1
                    tcal_i = tcal(i);
                else
                    tcal_i = tcal;
                end
                ntcal_i = 1;
                trans_i = trans(:,i);
                ntrans_i = 1;
                if ntransirf > 1
                    irf_i = irf(:,i);
                else
                    irf_i = irf;
                end
                ntransirf_i = 1;
                if ntransbck > 1
                    bck_i = bck(:,i);
                else
                    bck_i = bck;
                end
                ntransbck_i = 1;
                if nfitstart > 1
                    fit_start_i = fit_start(i);
                else
                    fit_start_i = fit_start;
                end
                nfitstart_i = 1;
                if nfitend > 1
                    fit_end_i = fit_end(i);
                else
                    fit_end_i = fit_end;
                end
                nfitend_i = 1;
                if ntranswgt > 1
                    wgts_i = wgts(:,i);
                else
                    wgts_i = wgts;
                end
                if npsi > 1
                    psi_i = psi(i);
                else
                    psi_i = psi;
                end
                npsi_i = 1;
                if ntaur > 1
                    taur_i = taur(i);
                else
                    taur_i = taur;
                end
                ntaur_i = 1;
                if nargout > 1
                    [results(i),spa_results(i)] = ...
                        fit_tcspc_dks_nlls_subroutine...
                        (tcal_i,ntcal_i,trans_i,ndata,ntrans_i,...
                        irf_i,ntransirf_i,bck_i,ntransbck_i,fit_start_i,...
                        nfitstart_i,fit_end_i,nfitend_i,param,nparam,...
                        ntau,ntheta,fxparam,gparam,spa_indx,...
                        spa_low,spa_high,spa_size,lb,nlb,ub,nub,wgts_i,wgtflag,...
                        psi_i,npsi_i,taur_i,ntaur_i,options,convflag,...
                        itpType,parflag,"Simultaneous");
                else
                    [results(i)] = ...
                        fit_tcspc_dks_nlls_subroutine...
                        (tcal_i,ntcal_i,trans_i,ndata,ntrans_i,...
                        irf_i,ntransirf_i,bck_i,ntransbck_i,fit_start_i,...
                        nfitstart_i,fit_end_i,nfitend_i,param,nparam,...
                        ntau,ntheta,fxparam,gparam,spa_indx,...
                        spa_low,spa_high,spa_size,lb,nlb,ub,nub,wgts_i,wgtflag,...
                        psi_i,npsi_i,taur_i,ntaur_i,options,convflag,...
                        itpType,parflag,"Simultaneous");
                end
                
            end
        end
    else
        % Generic workflow
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
                spa_size,"UniformOutput",false);
            nspadims = numel(spa_low);
            [spa_coord{1:nspadims}] = ndgrid(spagrids{:});
            %% Initialize SPA results
            sparesultsgrid = cellfun(@(x) 1:x,spa_size,"UniformOutput",false);
            spa_results_l(sparesultsgrid{:}) = results;
            nevals = numel(spa_results_l);
            % Fix all parameter values in SPA range
            fxparam([spa_indx{:}]) = true;
            if parflag
                parfor i = 1:nevals
                    % Return parameters to fix
                    spaparamvals = cellfun(@(x) x(i),spa_coord);
                    % Assign to param & fix
                    param_l = param;
                    spa_indx_l = spa_indx;
                    param_l([spa_indx_l{:}]) = spaparamvals;
                    % Recursively call this routine
                    spa_results_l(i) = fit_tcspc_dks_nlls_subroutine...
                        (tcal,ntcal,trans,ndata,ntrans,...
                        irf,ntransirf,bck,ntransbck,fit_start,nfitstart,fit_end,...
                        nfitend,param_l,nparam,ntau,ntheta,fxparam,gparam,spa_indx,...
                        spa_low,spa_high,spa_size,lb,nlb,ub,nub,wgts,wgtflag,...
                        psi,npsi,taur,ntaur,options,convflag,itpType,parflag,...
                        typeflag);
                end
                
            else
                for i = 1:nevals
                    % Return parameters to fix
                    spaparamvals = cellfun(@(x) x(i),spa_coord);
                    % Assign to param & fix
                    param([spa_indx{:}]) = spaparamvals;
                    % Recursively call this routine
                    spa_results_l(i) = fit_tcspc_dks_nlls_subroutine...
                        (tcal,ntcal,trans,ndata,ntrans,...
                        irf,ntransirf,bck,ntransbck,fit_start,nfitstart,fit_end,...
                        nfitend,param,nparam,ntau,ntheta,fxparam,gparam,spa_indx,...
                        spa_low,spa_high,spa_size,lb,nlb,ub,nub,wgts,wgtflag,...
                        psi,npsi,taur,ntaur,options,convflag,itpType,parflag,typeflag);
                end
            end
            spa_results.spagrid = spa_results_l;
        end
    end

%     %% Prepare engine arguments - Function handle
%     % DATA - Structure that contains all relevant parameters for fit
%     %%% TCAL - Bin width of data
%     %%% NTCAL - Number of bin widths
%     %%% TRANS - Data to compare to model
%     %%% NDATA - Number of bins in data
%     %%% NTRANS - Number of curves in model
%     %%% IRF - IRF data to use in convolution
%     %%% ITPIRF - Interpolated IRFs for color shifting
%     itpirf = interp_irf(tcal,ntcal,irf,ndata,ntransirf,itpType);
%     %%% NTRANSIRF - Number of IRFs
%     %%% BCK - Background data to use for artifact control
%     %%% NTRANSBCK - Number of backgrounds
%     %%% WGT - Square root of weights to use when provided
%     switch wgtflag
%         case "Gaussian - Data"
%             sqwgts = sqrt(gauss_wgts_data(trans));
%         case "Custom"
%             sqwgts = sqrt(wgts);
%         otherwise
%             sqwgts = wgts;
%     end
%     %%% NDATAWGT - Number of data points in weights
%     %%% NTRANSWGT - Number of curves in weights
%     [ndatawgt,ntranswgt] = size(sqwgts);
%     %%% WGTFLAG - Flag for dynamic computation of weights
%     %%% FIT_START - Where to start fitting
%     %%% NFITSTART - Number of fit start values
%     %%% FIT_END - Where to end fitting
%     %%% NFITEND - Number of fit end values
%     %%% PMAIN - Parameter vector for entire kernel
%     %%% MPOINT - Pointer matrix of values for aligning PMAIN
%     [pmain,mpoint,lbmain,ubmain] = get_pmain_mpoint_nlls...
%         (param,nparam,fxparam,gparam,lb,nlb,ub,nub,ntrans);
%     %%% FREEIDX - Free, unique indices in mpoint for indexing
%     freeidx = unique(mpoint(mpoint > 0));
%     %%% NPARAM - Number of parameters in model
%     %%% NTAU - Number of lifetimes in model
%     %%% NTHETA - Number of anisotropies in model
%     %%% FXPARAM - Fixed parameters
%     %%% GPARAM - Global parameters
%     %%% TAUR - Reference lifetimes to use for standard based convolution
%     %%% NTAUR - Number of reference lifetimes
%     %%% PSI - Polarization angle at which data were acquired
%     %%% NPSI - Number of polarization angles at which data were acquired
%     %%% CONVFLAG - Flag indicating which convolution routine to use
%     data = struct("tcal",tcal,"ntcal",ntcal,"trans",trans,"ndata",ndata,...
%         "ntrans",ntrans,"irf",irf,"itpirf",itpirf,"ntransirf",ntransirf,...
%         "bck",bck,"ntransbck",ntransbck,"wgt",sqwgts,"ndatawgt",ndatawgt,...
%         "ntranswgt",ntranswgt,"wgtflag",wgtflag,"fit_start",fit_start,...
%         "nfitstart",nfitstart,"fit_end",fit_end,"nfitend",nfitend,...
%         "pmain",pmain,"mpoint",mpoint,"freeidx",freeidx,"nparam",nparam,...
%         "ntau",ntau,"ntheta",ntheta,"fxparam",fxparam,"gparam",gparam,...
%         "taur",taur,"ntaur",ntaur,"psi",psi,"npsi",npsi,"convflag",convflag);
%     %%% F - Weighted residuals buffer for gated curves
%     [F,J] = prep_engine_buffers(fit_start,nfitstart,fit_end,nfitend,...
%         ntrans,numel(pmain));
%     fun = @(x) fit_tcspc_dks_lsqnonlin_objfun(x,data,F,J);
%     %% Prepare engine arguments - Other variables
%     x0 = pmain(freeidx);
%     if nlb > 0
%         lb_engine = lbmain(freeidx);
%     else
%         lb_engine = lb;
%     end
% 
%     if nub > 0
%         ub_engine = ubmain(freeidx);
%     else
%         ub_engine = ub;
%     end
%     %% Run engine for main routine
%     [x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin...
%         (fun,x0,lb_engine,ub_engine,options);
%     %% Organize outputs 
%     results = fit_tcspc_dks_lsqnonlin_output(x,resnorm,residual,exitflag,...
%         output,lambda,jacobian,data);
%     spaflag = ~isempty(spa_indx);
%     if spaflag && nargout > 1
%         %% Generate SPA grid
%         spagrids = cellfun(@(x,y,z) linspace(x,y,z),spa_low,spa_high,...
%             spa_size,"UniformOutput",false);
%         nspadims = numel(spa_low);
%         [spa_coord{1:nspadims}] = ndgrid(spagrids{:});
%         %% Initialize SPA results
%         sparesultsgrid = cellfun(@(x) 1:x,spa_size,"UniformOutput",false);
%         spa_results(sparesultsgrid{:}) = results;
%         nevals = numel(spa_results);
%         % Fix all parameter values in SPA range
%         fxparam([spa_indx{:}]) = true;
%         for i = 1:nevals
%             % Return parameters to fix
%             spaparamvals = cellfun(@(x) x(i),spa_coord);
%             % Assign to param & fix
%             param([spa_indx{:}]) = spaparamvals;
%             % Recursively call this routine
%             spa_results(i) = fit_tcspc_dks_nlls_subroutine...
%                 (tcal,ntcal,trans,ndata,ntrans,...
%                 irf,ntransirf,bck,ntransbck,fit_start,nfitstart,fit_end,...
%                 nfitend,param,nparam,ntau,ntheta,fxparam,gparam,spa_indx,...
%                 spa_low,spa_high,spa_size,lb,nlb,ub,nub,wgts,wgtflag,...
%                 psi,npsi,taur,ntaur,options,convflag,itpType);
%         end
%     end
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