function trans = sim_tcspc_dks(tcal,irf,bck,param,ntau,ncurves,args)
    %SIM_TCSPC_DKS Simulates decays acquired using TCSPC apparatus
    arguments
        tcal (1,1) {mustBePositive,mustBeFinite}
        irf (:,1) {mustBeNonnegative,mustBeInteger}
        bck (:,1) {mustBeNonnegative,mustBeInteger}
        param (:,1) {mustBeNumeric}
        ntau (1,1) {mustBeNonnegative,mustBeInteger}
        ncurves (1,1) {mustBeNonnegative,mustBeInteger} = 10000
        args.SimulationMethod {mustBeMember(args.SimulationMethod,...
            "Photons")} = "Photons"
        args.ConvolutionMethod {mustBeMember(args.ConvolutionMethod,...
            "FFT")} = "FFT"
        args.InterpolationOption {mustBeMember(args.InterpolationOption,...
            ["smoothingspline","cubicinterp"])} = "smoothingspline"
    end
    %% Check that number of irf/background curves match
    nirfbins = numel(irf);
    nbckbins = numel(bck);
    assert(nirfbins == nbckbins);
    %% Check for parameter size agreement
    nparam = numel(param);
    ntaucomp = 2 * ntau;
    nexpparam = ntaucomp + 4;
    assert(nparam == nexpparam);
    switch args.SimulationMethod
        case "Photons"
            trans = sim_tcspc_dks_photons(tcal,irf,bck,param,ntau,...
                ncurves,args.ConvolutionMethod,args.InterpolationOption);
        otherwise
            error("Unsupported simulation approach")
    end
end

function trans = sim_tcspc_dks_photons(tcal,irf,bck,param,ntau,ncurves,cMethod,iOpt)
    % Subroutine for photon simulation workflow
    %% Housekeeping
    nbins = numel(irf);
    T = tcal * nbins;
    %% Check that all amplitude terms are in photon counts
    ntaucomp = 2 * ntau;
    ntaucounts = param(1:2:ntaucomp);
    nlincounts = vertcat(ntaucounts,param(end - 3:end));
    mustBeInteger(nlincounts);
    mustBeNonnegative(nlincounts);
    %% Simulate exponential photons
    expphots_cell = cell(ntau,1);
    taus = param(2:2:ntaucomp);
    expphots_cell = fill_expphots_cell(expphots_cell,taus,ntau,ntaucounts,ncurves);
    ntotaltaucounts = sum(ntaucounts);
    expphots = zeros(ntotaltaucounts,ncurves,class(expphots_cell{1}));
    expphots = fill_expphots(expphots,expphots_cell,ntaucounts,ntau);
    %expcounts = cell2mat(expphots_cell);
    %% Shift IRF
    Q = param(end - 3);
    irf = shift_irf(irf,tcal,Q,iOpt);
    %% Normalize IRF
    irf = irf / sum(irf);
    %% Sample from shifted IRF using ITS
    nscphots = param(end - 2);
    scphots = inverse_transform_sampling(tcal,irf,iOpt,nscphots,ncurves);
    %% Sample from background using ITS
    nbckphots = param(end - 1);
    bckphots = inverse_transform_sampling(tcal,bck,iOpt,nbckphots,ncurves);
    %% Sample from uniform background
    noffsphots = param(end);
    offsphots = random("Uniform",0,T,noffsphots,ncurves);
    %% Put photons into big matrix
    photmatrix = vertcat(expphots,scphots,bckphots,offsphots);
    %% Bin photons into curves
    trans = bin_photons(tcal,nbins,ncurves,photmatrix);
    %% Convolve photons with PDF of IRF
    switch cMethod
        case "FFT"
            trans = round(abs(ifft(fftshift(fft(trans,[],1) .* fft(irf,[],1)),[],1)));
        otherwise
            error("Unsupported convolution!")
    end
    %% Apply Poisson weighting
    trans = random("Poisson",trans);
end

function expphots_cell = fill_expphots_cell(expphots_cell,taus,ntau,ntaucounts,ncurves)
    % Fills cell array with exponential counts
    for i = 1:ntau
        expphots_cell{i} = random("Exponential",taus(i),ntaucounts(i),ncurves);
    end
end

function expphots = fill_expphots(expphots,expphots_cell,ntaucounts,ntau)
    offs = 0;
    for i = 1:ntau
        localcounts = ntaucounts(i);
        expphots((offs + 1):(offs + localcounts),:) = expphots_cell{i};
        offs = offs + localcounts;
    end
end

