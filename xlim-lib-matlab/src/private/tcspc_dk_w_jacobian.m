function J = tcspc_dk_w_jacobian(tcal,ndata,irf,bck,wgt,wgtflag,param,...
    nparam,ntau,ntheta,taur,psi,convflag)
    %TCSPC_DK_W_JACOBIAN Computes weighted Jacobian from given curve

    %% Compute model
    if ntheta > 0
        
    else
        [~,Jalpha,Jtau] = tcspc_dk_flim_data_model...
            (tcal,ndata,irf,param,ntau,taur,convflag);
    end
    %% Account for artifacts
    Jsc = irf; Jbck = bck; Jz = ones(size(irf));
    J = horzcat(Jalpha,Jtau,Jsc,Jbck,Jz);
    J(:,1:2:(2 * ntau)) = Jalpha;
    J(:,2:2:(2 * ntau)) = Jtau;
    %% Determine which weights to use
    switch wgtflag
        otherwise
            sqwgts = wgt;
    end
    J = sqwgts .* J;
end

