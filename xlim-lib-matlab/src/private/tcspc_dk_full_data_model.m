function M = tcspc_dk_full_data_model(tcal,ndata,irf,bck,param,nparam,ntau,ntheta,taur,psi,convflag)
    %TCSPC_DK_FULL_DATA_MODEL Computes full data model for inputs
    %% Compute pure model
    if ntheta > 0
    else
        I = tcspc_dk_flim_data_model(tcal,ndata,irf,param,ntau,taur,convflag);
    end
    %% Account for artifacts
    S = param(nparam - 2); V = param(nparam - 1); Z = param(nparam);
    M = I + (S * irf) + (V * bck) + Z;
end

