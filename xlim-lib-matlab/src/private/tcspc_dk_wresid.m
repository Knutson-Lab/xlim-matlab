function wresid = tcspc_dk_wresid(tcal,trans,ndata,irf,bck,wgt,wgtflag,param,...
    nparam,ntau,ntheta,taur,psi,convflag)
    %TCSPC_DK_WRESID Computes weighted residual for single curve
    %% Compute pure model
    if ntheta > 0
    else
        I = tcspc_dk_flim_data_model(tcal,ndata,irf,param,ntau,taur,convflag);
    end
    %% Account for artifacts
    S = param(nparam - 2); V = param(nparam - 1); Z = param(nparam);
    M = I + (S * irf) + (V * bck) + Z;
    %% Determine weights to use
    switch wgtflag
        otherwise
            sqwgts = wgt;
    end
    wresid = sqwgts .* (M - trans);
end

