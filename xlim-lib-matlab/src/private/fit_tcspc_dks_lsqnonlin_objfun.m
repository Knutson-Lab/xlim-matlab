function [F,J] = fit_tcspc_dks_lsqnonlin_objfun(x0,data,F,J_full)
    %FIT_TCSPC_DKS_NLLS_OBJFUN - Objective function for NLLS routine driven
    %via lsqnonlin
    %% Update pmain with parameters
    data.pmain(data.freeidx) = x0;
    %% Fill buffer with data
    F = fill_F(F,data);
    if nargout > 1
        J_full = fill_J_full(J_full,data);
        J = J_full(:,data.freeidx);
    end
end

function F = fill_F(F,data)
    %% Initialize starting gate for buffer
    sc = 0;
    for i = 1:data.ntrans
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
        % Determine weight
        if data.ntranswgt > 1
            wgt_i = data.wgt(:,i);
        else
            wgt_i = data.wgt;
        end
        % Determine fit start/stop for current curve
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
        % Calculate weighted residual for given curve
        wresid = tcspc_dk_wresid(tcal_i,data.trans(:,i),data.ndata,...
            irf_i,bck_i,wgt_i,data.wgtflag,param,data.nparam,data.ntau,...
            data.ntheta,taur_i,psi_i,data.convflag);
        % Put weighted residual back into buffer
        ec = fit_end_i - fit_start_i + 1;
        F((sc + 1):(sc + ec)) = wresid(fit_start_i:fit_end_i);
        sc = sc + ec;
    end
    
end

function J_full = fill_J_full(J_full,data)
    sc = 0;
    for i = 1:data.ntrans
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
        % Determine weight
        if data.ntranswgt > 1
            wgt_i = data.wgt(:,i);
        else
            wgt_i = data.wgt;
        end
        % Determine fit start/stop for current curve
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
        % Calculate weighted Jacobian for given curve
        wJ = tcspc_dk_w_jacobian(tcal_i,data.ndata,irf_i,bck_i,wgt_i,...
            data.wgtflag,param,data.nparam,data.ntau,data.ntheta,taur_i,...
            psi_i,data.convflag);
        % Put weighted residual back into buffer
        ec = fit_end_i - fit_start_i + 1;
        J_full((sc + 1):(sc + ec),mpoint_i(freepar)) = ...
            J_full((sc + 1):(sc + ec),mpoint_i(freepar)) + ...
            wJ(fit_start_i:fit_end_i,freepar);
        sc = sc + ec;
    end
end

