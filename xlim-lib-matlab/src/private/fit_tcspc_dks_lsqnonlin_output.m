function results = fit_tcspc_dks_lsqnonlin_output(x,resnorm,residual,...
    exitflag,output,lambda,jacobian,data)
    %FIT_TCSPC_DKS_LSQNONLIN_OUTPUT Organizes results from engine into more
    %legible output

    %% Put most recent fits back into pmain
    data.pmain(data.freeidx) = x;
    %% Organize results - Curves
    curveresults(1:data.ntrans,1) = struct("fitted",zeros(data.ndata,1),...
        "residuals",zeros(data.ndata,1),"param",zeros(data.nparam,1),...
        "sterr",zeros(data.nparam,1),"X2",0,"df",0);
    curveresults = fill_curve_results(curveresults,data,residual,jacobian);
    results.curves = curveresults;
    %% Organize results - Engine
    results.X2 = resnorm;
    results.df = numel(residual) - numel(x);
    results.residuals = residual;
    results.exitflag = exitflag;
    results.output = output;
    results.lambda = lambda;
    results.jacobian = jacobian;
end

function curveresults = fill_curve_results(curveresults,data,residual,jacobian)
    sc = 0;
    J(:,data.freeidx) = jacobian;
    for i = 1:data.ntrans
        %% Extract necessary variables
        % Determine time axis 
        if data.ntcal > 1
            tcal_i = data.tcal(i);
        else
            tcal_i = data.tcal;
        end
        tax_i = (tcal_i / 2):tcal_i:(data.ndata * tcal_i);
        % Determine IRF
        if data.ntransirf > 1
            irf_i = data.irf(:,i);
        else
            irf_i = data.irf;
        end
        % Determine background 
        if data.ntransbck > 1
            bck_i = data.bck(:,i);
        else
            bck_i = data.bck;
        end
        % Get channels
        if data.nfitstart > 1
            fit_start_i = data.fit_start(i);
        else
            fit_start_i = data.fit_start;
        end
        if data.nfitend > 1
            fit_end_i = data.fit_end(i);
        else
            fit_end_i = data.fit_end;
        end
        % Determine parameter
        mpoint_i = data.mpoint(:,i);
        param = data.pmain(abs(mpoint_i));
        freepar = ~data.fxparam;
        % Determine reference lifetime
        if data.ntaur > 1
            taur_i = data.taur(i);
        else
            taur_i = data.taur;
        end
        % Determine polarization
        if data.npsi > 1
            psi_i = data.psi(i);
        else
            psi_i = data.psi;
        end
        % Shift IRF by Q
        q = param(data.nparam - 3);
        if data.ntransirf > 1
            itpirf_i = data.itpirf{i};
        else
            itpirf_i = data.itpirf;
        end
        irf_i(:) = abs(round(itpirf_i(tax_i + (q * tcal_i))));
        %% Compute fitted curve
        curveresults(i).fitted = tcspc_dk_full_data_model...
            (tcal_i,data.ndata,irf_i,bck_i,param,data.nparam,data.ntau,...
            data.ntheta,taur_i,psi_i,data.convflag);
        %% Extract residuals
        ec = fit_end_i - fit_start_i + 1;
        curveresults(i).residuals(fit_start_i:fit_end_i) = ...
            residual((sc + 1):(sc + ec));
        curveresults(i).param = param;
        %% Compute standard error around parameter
        wJ = J((sc + 1):(sc + ec),mpoint_i(freepar));
        C = (wJ' * wJ) ^ (-1);
        %% Compute computed variance of residuals
        nfreepar = nnz(freepar);
        curveresults(i).df = (fit_end_i - fit_start_i - nfreepar);
        curveresults(i).X2 = curveresults(i).residuals' * ...
            curveresults(i).residuals;
        curveresults(i).covar = C * (curveresults(i).X2 / curveresults(i).df);
        curveresults(i).sterr = sqrt(diag(curveresults(i).covar));
        sc = sc + ec;
    end
end