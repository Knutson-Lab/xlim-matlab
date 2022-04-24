function results = fit_tcspc_dks_nlls(tcal,trans,irf,bck,fit_start,fit_end,...
    param,fxparam,gparam,ntau,ntheta,wgts,lb,ub,...
    spa_indx,spa_low,spa_high,spa_size,args)
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
        wgts (:,:) {mustBeNonnegative,mustBeFinite} = gauss_wgts(trans)
        lb (:,1) {mustBeNumeric} = []
        ub (:,1) {mustBeNumeric} = []
    end

    arguments (Repeating)
        spa_indx (1,1) {mustBePositive,mustBeInteger}
        spa_low (1,1) {mustBeNumeric,mustBeFinite}
        spa_high (1,1) {mustBeNumeric,mustBeFinite}
        spa_size (1,1) {mustBePositive,mustBeInteger}
    end

    arguments
        args.PolarizationAngle (:,1) {mustBeInRange...
            (args.PolarizationAngle,0,90)}
        args.ReferenceLifetime (:,1) {mustBePositive,mustBeFinite}
        args.EngineOptions (1,1) optim.options.Lsqnonlin = optimoptions("lsqnonlin")
        args.ConvolutionMethod {mustBeMember(args.ConvolutionMethod,...
            "FFT")} = "FFT"
        args.InterpolationOption {mustBeMember(args.InterpolationOption,...
            ["smoothingspline","cubicinterp"])} = "smoothingspline"
    end

    %% Valdiate required arguments
    validate_req_args(tcal,trans,irf,bck,fit_start,fit_end,param,...
        fxparam,gparam,ntau,ntheta)
    %% Validate optional arguments
    validate_opt_args(trans,param,wgts,lb,ub)
    %% Validate repeating arguments
    validate_repeat_args(param,spa_indx,spa_low,spa_high)
    %% Validate N/V arguments
    validate_nv_args(tcal,irf,fxparam,gparam,args);
    if ~isfield(args,"ReferenceLifetime")
        args.ReferenceLifetime = [];
    end
    if ~isfield(args,"PolarizationAngle")
        args.PolarizationAngle = [];
    end
    %% For now, split based on convolution routine
    switch args.ConvolutionMethod
        case "FFT"
            results = fit_tcspc_dks_nlls_fft(tcal,trans,irf,bck,...
                fit_start,fit_end,param,fxparam,gparam,ntau,ntheta,wgts,...
                lb,ub,spa_indx,spa_low,spa_high,spa_size,...
                args.PolarizationAngle,args.ReferenceLifetime,...
                args.EngineOptions,args.InterpolationOption);
    end
end

function results = fit_tcspc_dks_nlls_fft(tcal,trans,irf,bck,fit_start,...
    fit_end,param,fxparam,gparam,ntau,ntheta,wgts,lb,ub,spa_indx,...
    spa_low,spa_high,spa_size,psi,taur,engineOpts,itpOpts)
    % Subroutine for NLLS fitting with a FFT convolution
    %% Prepare time axes
    t = prep_t_axis(tcal,trans);
    is_t_scalar = size(t,2) == 1;
    %% Prepare IRF interpolations
    irfitps = prep_irf_itps(tcal,irf,itpOpts);
    is_irf_scalar = isscalar(irfitps);
    nirfs = numel(irfitps);
    %% Prepare parameter variables
    [ndata,ncurves] = size(trans);
    [pmain,mpoint,lbmain,ubmain] = get_pmain_mpoint_nlls...
        (param,fxparam,gparam,lb,ub,ncurves);
    abs_mpoint = abs(mpoint);
    freeindx = unique(mpoint(mpoint > 0));
    x0 = pmain(freeindx);
    if ~isempty(lbmain)
        lb_eng = lbmain(freeindx);
    else
        lb_eng = [];
    end
    if ~isempty(ubmain)
        ub_eng = ubmain(freeindx);
    else
        ub_eng = [];
    end
    
    %% Prepare specific indexing variables
    B_alpha = zeros(1,ncurves,ntau);
    alpha_idx = 1:2:(2 * ntau);
    B_tau = zeros(1,ncurves,ntau);
    tau_idx = 2:2:(2 * ntau);
    theta_offs = (2 * ntau) + 1;
    B_beta = zeros(1,ncurves,1,ntheta);
    beta_idx = theta_offs:2:(theta_offs + (2 * ntheta));
    B_theta = zeros(1,ncurves,1,ntheta);
    theta_idx = (theta_offs + 1):2:(theta_offs + (2 * ntheta));
    beta_flag = ntheta > 0;
    nparam = numel(param);
    q_idx = nparam - 3;
    sc_idx = nparam - 2;
    bck_idx = nparam - 1;
    offs_idx = nparam;
    %% Prepare engine buffers
    nfreepar = numel(x0);
    [B_resid,B_jac] = prep_engine_buffers(fit_start,fit_end,ncurves,nfreepar);
    is_fit_start_scalar = isscalar(fit_start);
    is_fit_end_scalar = isscalar(fit_end);
    B_jac_mod = zeros(ndata,ncurves,nparam);
    %% Take square root of weights
    sqwgts = prep_wgts(fit_start,fit_end,ncurves,wgts);
    %sqwgts = sqrt(wgts);
    %% Put reference lifetime in appropriate container
    taur_flag = ~isempty(taur);
    if taur_flag
        B_taur = zeros(1,ncurves,1);
        B_taur(1,:,1) = taur;
    else
        B_taur = [];
    end 
    %% Precompute eta
    eta = (3 * (cos(psi) .^ 2)) - 1;
    %% Specify objective function 
    freepar = ~fxparam;
    fun = @(x) nlls_fft_obj(x,tcal,t,is_t_scalar,trans,ncurves,...
        irf,nirfs,is_irf_scalar,irfitps,bck,fit_start,is_fit_start_scalar,...
        fit_end,is_fit_end_scalar,...
        pmain,freepar,mpoint,abs_mpoint,freeindx,B_alpha,alpha_idx,...
        B_tau,tau_idx,B_beta,beta_idx,B_theta,theta_idx,beta_flag,...
        q_idx,sc_idx,bck_idx,offs_idx,B_resid,B_jac,B_jac_mod,sqwgts,...
        B_taur,taur_flag,eta);
    [x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin...
        (fun,x0,lb_eng,ub_eng,engineOpts);
    x;
end

function sqwgts = prep_wgts(fit_start,fit_end,ncurves,wgts)
    chanwindow = fit_end - fit_start + 1;
    if isscalar(fit_start) && isscalar(fit_end)
        ntotalbins = ncurves * chanwindow;
    else
        ntotalbins = sum(chanwindow);
    end
    B_wgts = zeros(ntotalbins,1);
    B_wgts = fill_F(B_wgts,ncurves,fit_start,isscalar(fit_start),fit_end,...
        isscalar(fit_end),wgts);
    sqwgts = sparse(ntotalbins,ntotalbins);
    sqwgts(logical(eye(ntotalbins))) = B_wgts;
    sqwgts = sqrtm(sqwgts);
end

function [F,J] = nlls_fft_obj(x0,tcal,t,is_t_scalar,trans,ntrans,...
    irf,nirf,is_irf_scalar,irfitps,bck,fit_start,is_fit_start_scalar,...
    fit_end,is_fit_end_scalar,pmain,freepar,mpoint,abs_mpoint,...
    freeidx,B_alpha,alpha_idx,B_tau,tau_idx,B_beta,beta_idx,B_theta,...
    theta_idx,beta_flag,q_idx,sc_idx,bck_idx,offs_idx,F,J,B_jac_mod,...
    sqwgts,B_taur,taur_flag,eta)
    %% Update pmain with guesses
    pmain(freeidx) = x0;
    %% Update IRF with new shift value
    q = pmain(abs_mpoint(q_idx,:)) .* tcal;
    irf = fill_shift_irf(irf,irfitps,is_irf_scalar,nirf,t,is_t_scalar,q);
    %% Calculate model values
    mod = model_matrix(t,irf,bck,pmain,abs_mpoint,B_alpha,alpha_idx,...
        B_tau,tau_idx,B_beta,beta_idx,B_theta,theta_idx,beta_flag,...
        sc_idx,bck_idx,offs_idx,B_taur,taur_flag,eta);
    %% Compute weighted residuals
    w_resid = mod - trans;
    %w_resid = sqwgts .* (mod - trans);
    %% Fill buffer
    F = fill_F(F,ntrans,fit_start,is_fit_start_scalar,fit_end,...
        is_fit_end_scalar,w_resid);
    F = sqwgts * F;
    if nargout > 1
        % Calculate jacobian
        B_jac_mod = model_jacobian(B_jac_mod,t,irf,bck,pmain,abs_mpoint,...
            B_alpha,alpha_idx,B_tau,tau_idx,B_beta,beta_idx,B_theta,...
            theta_idx,beta_flag,sc_idx,bck_idx,offs_idx,B_taur,taur_flag,...
            eta);
        % Compute weighted jacobian
        %B_jac_mod = sqwgts .* B_jac_mod;
        % Fill buffer
        J = fill_J(J,ntrans,fit_start,is_fit_start_scalar,fit_end,...
            is_fit_end_scalar,mpoint,freepar,B_jac_mod);
        J = sqwgts * J;
    end
end

function J = fill_J(J,ntrans,fit_start,is_fit_start_scalar,fit_end,...
    is_fit_end_scalar,mpoint,freepar,B_jac_mod)
    sc_g = 0;
    for i = 1:ntrans
        if is_fit_start_scalar
            sc_l = fit_start;
        else
            sc_l = fit_start(i);
        end
        if is_fit_end_scalar
            ec_l = fit_end;
        else
            ec_l = fit_end(i);
        end
        win_l = ec_l - sc_l + 1;
        indx_loc = mpoint(freepar,i);
        J((sc_g + 1):(sc_g + win_l),indx_loc) = B_jac_mod(sc_l:ec_l,i,freepar);
        sc_g = sc_g + win_l;
    end
end

function F = fill_F(F,ntrans,fit_start,is_fit_start_scalar,...
    fit_end,is_fit_end_scalar,w_resid)
    sc_g = 0;
    for i = 1:ntrans
        if is_fit_start_scalar
            sc_l = fit_start;
        else
            sc_l = fit_start(i);
        end
        if is_fit_end_scalar
            ec_l = fit_end;
        else
            ec_l = fit_end(i);
        end
        win_l = ec_l - sc_l + 1;
        F((sc_g + 1):(sc_g + win_l)) = w_resid(sc_l:ec_l);
        sc_g = sc_g + win_l;
    end
    
end

function J_mod = model_jacobian(J_mod,t,irf,bck,pmain,abs_mpoint,B_alpha,...
    alpha_idx,B_tau,tau_idx,B_beta,beta_idx,B_theta,theta_idx,...
    beta_flag,sc_idx,bck_idx,offs_idx,taur,taur_flag,eta)
    %% Prepare non-artifact terms
    B_alpha(1,:,:) = pmain(abs_mpoint(alpha_idx,:));
    B_tau(1,:,:) = pmain(abs_mpoint(tau_idx,:));
    B_beta(1,:,:) = pmain(abs_mpoint(beta_idx,:));
    B_theta(1,:,:) = pmain(abs_mpoint(theta_idx,:));
    trans_irf = fft(irf,[],1);
    %% Compute initial derivative
    D1tmp = D(t,B_tau);
    Dalpha1 = D_tilde(t,B_tau,trans_irf);
    Dtau1 = abs(ifft(fftshift(fft(B_alpha .* (t ./ (B_tau .^ 2)) .* D1tmp,[],1) .* trans_irf),[],1));
    if taur_flag
        Dalpha1 = (Dalpha1 ./ taur) - (Dalpha1 ./ B_tau);
        Dtau1 = (Dtau1 ./ taur) + ...
            abs(ifft(fftshift(fft((B_alpha .* D1tmp .* ((B_tau - t) ./ (B_tau .^ 3))),[],1)),[],1));
    end
    B_tau_eff = anisotropy_vals(B_tau,B_theta);
    D2tmp = D(t,B_tau_eff); 
    Dalpha2 = D_tilde(t,B_tau_eff,trans_irf);
    Dtau2 = abs(ifft(fftshift(fft(B_alpha .* (t ./ (B_tau .^ 2)) .* D2tmp,[],1) .* trans_irf),[],1));
    Dbeta = Dalpha2;
    Dtheta = abs(ifft(fftshift(fft(B_beta .* (t ./ (B_theta .^ 2)) .* D2tmp,[],1) .* trans_irf),[],1));
    if taur_flag
        Dalpha2 = (Dalpha2 ./ taur) - (Dalpha2 ./ B_tau);
        Dtau2 = (Dtau2 ./ taur) + ...
            abs(ifft(fftshift(fft(B_alpha .* (D2tmp .* ((B_tau - t) ./ (B_tau .^ 3))),[],1) .* trans_irf),[],1));
        Dbeta = (Dbeta ./ taur) - (Dbeta ./ B_theta);
        Dtheta = (Dtheta ./ taur) + ...
            abs(ifft(fftshift(fft(B_beta .* (D2tmp .* ((B_theta - t) ./ (B_theta .^ 3))),[],1) .* trans_irf),[],1));
    end
    if beta_flag
        J_mod(:,:,alpha_idx) = (1 / 3) * Dalpha1 + ...
            (eta / 3) .* sum(B_beta .* Dalpha2, 4);
        J_mod(:,:,tau_idx) = (1 / 3) * Dtau1 + ...
            (eta / 3) .* sum(B_beta .* Dtau2, 4);
        J_mod(:,:,beta_idx) = permute...
            ((eta / 3) .* sum(B_alpha .* Dbeta,3),[1 2 4 3]);
        J_mod(:,:,theta_idx) = permute...
            ((eta /3) .* sum(B_alpha .* Dtheta,3),[1 2 4 3]);
    else
        J_mod(:,:,alpha_idx) = Dalpha1;
        J_mod(:,:,tau_idx) = Dtau1;
        J_mod(:,:,beta_idx) = [];
        J_mod(:,:,theta_idx) = [];
    end
    %% Compute artifact derivatives
    J_mod(:,:,sc_idx) = irf;
    J_mod(:,:,bck_idx) = bck;
    J_mod(:,:,offs_idx) = 1;
end

function mat = model_matrix(t,irf,bck,pmain,abs_mpoint,B_alpha,...
    alpha_idx,B_tau,tau_idx,B_beta,beta_idx,...
    B_theta,theta_idx,beta_flag,sc_idx,bck_idx,offs_idx,taur,taur_flag,eta)
    %% Prepare non-artifact model
    B_alpha(1,:,:) = pmain(abs_mpoint(alpha_idx,:));
    B_tau(1,:,:) = pmain(abs_mpoint(tau_idx,:));
    B_beta(1,:,:) = pmain(abs_mpoint(beta_idx,:));
    B_theta(1,:,:) = pmain(abs_mpoint(theta_idx,:));
    trans_irf = fft(irf,[],1);
    D1 = D_tilde(t,B_tau,trans_irf);
    if taur_flag
        D1 = D_tilde_r(D1,B_tau,irf,taur);
    end
    B_tau_eff = anisotropy_vals(B_tau,B_theta);
    D2 = D_tilde(t,B_tau_eff,trans_irf);
    if taur_flag
        D2 = D_tilde_r(D2,B_tau_eff,irf,taur);
    end
    if beta_flag
        % Do anisotropy
        mat = I_tilde_psi(B_alpha,D1,eta,B_beta,D2);
    else
        % No anisotropy
        mat = I_tilde(B_alpha,D1);
    end
    %% Now account for artifacts
    S = pmain(abs_mpoint(sc_idx,:));
    V = pmain(abs_mpoint(bck_idx,:));
    Z = pmain(abs_mpoint(offs_idx,:));
    mat = M(mat,S,irf,V,bck,Z);
end

function irf = fill_shift_irf(irf,irfitps,irf_flag,nirf,t,t_flag,q)
    shifted_t = ((t') + q)';
    if irf_flag
        % Scalar IRF
        irf(:) = irfitps{1}(shifted_t(:,1));
    else
        % Multiple irfs
        for i = 1:nirf
            if t_flag
                irf(:,i) = irfitps{i}(shifted_t(:,1));
            else
                irf(:,i) = irfitps{i}(shifted_t(:,i));
            end
            
        end
    end
    irf = abs(round(irf));
end

function total_mat = M(I,S,irf,V,bck,Z)
    total_mat = I + S * irf + V * bck + Z;
end

function tau_theta_mat = anisotropy_vals(tau,theta)
    tau_theta_mat = ((1 ./ tau) + (1 ./ theta)) .^ (-1);
end

function conv_mat_tau = I_tilde_psi(alpha,Dtilde,eta,beta,Dtilde2)
    conv_mat_tau = (1 / 3) * sum(alpha .* Dtilde) + ...
        (eta / 3) * sum(beta .* sum(alpha .* Dtilde2,3),4);
end

function conv_mat_tau = I_tilde(amp,Dtilde)
    conv_mat_tau = sum(amp .* Dtilde,3);
end

function conv_mat_tau = D_tilde_r(Dtilde,tau,g,taur)
    conv_mat_tau = g + ((1 / taur) - (1 ./ tau)) .* Dtilde;
end

function conv_mat_tau = D_tilde(t,tau,trans_irf)
    conv_mat_tau = abs(ifft(fftshift(fft(D(t,tau),[],1) .* trans_irf),[],1));
end

function mat_tau = D(t,tau)
    mat_tau = exp(-t ./ tau);
end

function [B_resid,B_jac] = prep_engine_buffers(fit_start,fit_end,ncurves,nfreepar)
    chanwindow = fit_end - fit_start + 1;
    if isscalar(fit_start) && isscalar(fit_end)
        ntotalbins = ncurves * chanwindow;
    else
        ntotalbins = sum(chanwindow);
    end
    B_resid = zeros(ntotalbins,1);
    B_jac = zeros(ntotalbins,nfreepar);
end

function irfitps = prep_irf_itps(tcal,irf,itpOpts)
    ndata = size(irf,2);
    % Initialize cfit object by fitting first
    t = prep_t_axis(tcal,irf);
    is_t_scalar = size(t,2) == 1;
    ttemp = t(:,1);
    irfitps = repmat({fit(ttemp,irf(:,1),itpOpts)},ndata,1);
    for i = 1:ndata
        if is_t_scalar
            irfitps(i) = {fit(ttemp,irf(:,i),itpOpts)};
        else
            irfitps(i) = {fit(t(:,i),irf(:,i),itpOpts)};
        end
    end
end

function t = prep_t_axis(tcal,trans)
    nbins = size(trans,1);
    tstart = tcal / 2;
    tend = nbins * tcal;
    ntcal = numel(tcal);
    t = zeros(nbins,ntcal);
    for i = 1:ntcal
        t(:,i) = tstart(i):tcal(i):tend(i);
    end
end

function [pmain,mpoint,lbmain,ubmain] = ...
    get_pmain_mpoint_nlls(param,fxparam,gparam,lb,ub,ncurves)
    % Subroutine to compute parameter vector with pointer matrix
    %% Find mpoint first
    mpoint = get_mpoint(fxparam,gparam,ncurves);
    %% Return parameter vector
    Pm = repmat(param,1,ncurves);
    Mpa = abs(mpoint);
    mxind = max(Mpa,[],'all');
    pmain = zeros(mxind,1);
    pmain(Mpa) = Pm;
    if ~isempty(lb)
        Pm = repmat(lb,1,ncurves);
        lbmain = zeros(mxind,1);
        lbmain(Mpa) = Pm;
    else
        lbmain = [];
    end
    if ~isempty(ub)
        Pm = repmat(ub,1,ncurves);
        ubmain = zeros(mxind,1);
        ubmain(Mpa) = Pm;
    else
        ubmain = [];
    end
    
end

function validate_nv_args(tcal,irf,fxparam,gparam,args)
    % Validates that engine is compatible with config
    isqfx = fxparam(end - 3);
    anyglob = any(gparam);
    isobjgrad = args.EngineOptions.SpecifyObjectiveGradient;
    %% If there's a global parameter, objective gradient must be specified
    if anyglob && ~isobjgrad
        error("Global NLLS requires objective gradient")
    end
    %% If there's a free color shift and an objective gradient, cannot run
    if ~isqfx && isobjgrad
        error("Color shift estimate only works with numeric derivatives")
    end
    %% Check for polarization angle
    if isfield(args,"PolarizationAngle")
        nangles = numel(args.PolarizationAngle);
        assert(numel(tcal) == nangles || nangles == 1);
    end
    %% Check for reference lifetimes
    if isfield(args,"ReferenceLifetime")
        ntaus = numel(args.ReferenceLifetime);
        assert(size(irf,2) == ntaus || ntaus == 1)
    end
end

function validate_repeat_args(param,spa_indx,spa_low,spa_high)
    % Validates repeating arguments
    nparam = numel(param);
    assert(all(cellfun(@(x) x <= nparam,spa_indx)))
    assert(all(cellfun(@(x,y) x < y,spa_low,spa_high)));
end

function validate_opt_args(trans,param,wgts,lb,ub)
    % Validates sizes of optional positional arguments
    %% Housekeeping
    [nbins,ndata] = size(trans);
    %% Check weights
    [nwgtbins,nwgts] = size(wgts);
    assert((nwgtbins == nbins && (nwgts == ndata || ncurves == 1)) || isscalar(wgts))
    %% Check lower/upper bounds
    nparam = numel(param);
    validate_param_sz(lb,nparam);
    validate_param_sz(ub,nparam);
    assert(all(lb < ub))

end

function validate_req_args(tcal,trans,irf,bck,fit_start,fit_end,param,fxparam,gparam,ntau,ntheta)
    % Validates sizes of required arguments

    %% Establishing variables
    [nbins,ndata] = size(trans);
    %% Check sizes of data curves
    ntcal = numel(tcal);
    assert(ntcal == ndata || isscalar(tcal));
    validate_curve_sz(irf,nbins,ndata);
    validate_curve_sz(bck,nbins,ndata);
    %% Check dimensions of fit_start, fit_end
    nfitstart = numel(fit_start);
    nfitend = numel(fit_end);
    assert(nfitstart == nfitend || isscalar(fit_start) || isscalar(fit_end))
    assert(all(fit_start < fit_end) && all(fit_start <= nbins))
    %% Check sizes of parameters
    nparam = numel(param);
    nexpparam = (2 * ntau) + (2 * ntheta) + 4;
    assert(nparam == nexpparam);
    validate_param_sz(fxparam,nparam)
    validate_param_sz(gparam,nparam)
end

function validate_curve_sz(curvedata,nbins,ndata)
    [ncurvebins,ncurves] = size(curvedata);
    assert(ncurvebins == nbins && (ncurves == ndata || ncurves == 1))
end

function validate_param_sz(paramvec,nparam)
    noparam = numel(paramvec);
    if noparam > 0
        assert(nparam == noparam)
    end
    
end

function wgts = gauss_wgts(trans)
    % Subroutine that gets weights of 1 / N
    wgts = trans;
    nz = wgts > 0;
    wgts(nz) = 1 ./ wgts(nz);
    wgts(~nz) = 1;
end

