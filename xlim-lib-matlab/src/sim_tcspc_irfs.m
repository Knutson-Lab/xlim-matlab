function irfs = sim_tcspc_irfs(tcal,bck,param,ncurves,args)
    %SIM_TCSPC_IRFS Simulates instrument response functions from TCSPC
    %setup

    arguments
        tcal (1,1) {mustBePositive,mustBeFinite}
        bck (:,1) {mustBeNonnegative,mustBeInteger}
        param (:,1) {mustBeNumeric}
        ncurves (1,1) {mustBeNonnegative,mustBeInteger} = 10000
        args.SimulationMethod {mustBeMember(args.SimulationMethod,...
            "Photons")} = "Photons"
        args.InterpolationOption {mustBeMember(args.InterpolationOption,...
            ["smoothingspline","cubicinterp"])} = "smoothingspline"
    end

    %% Validate parameter size
    nparam = numel(param);
    assert(nparam == 5)
    %% Simulate IRF photons
    switch args.SimulationMethod
        case "Photons"
            irfs = sim_tcspc_irfs_photons(tcal,bck,param,ncurves,...
                args.InterpolationOption);
        otherwise
            error("Unsupported simulation method!")
    end
end

function irf = sim_tcspc_irfs_photons(tcal,bck,param,ncurves,iOpt)
    % Subroutine for simulating TCSPC IRFs using photon based approach
    %% Housekeeping
    nbins = numel(bck);
    T = tcal * nbins;
    %% Check that all amplitudes are nonnegative integers
    photparam = param([1,4:5]);
    mustBeNonnegative(photparam);
    mustBeInteger(photparam);
    %% Simulate Gaussian photons
    ngaussphots = param(1);
    mu = param(2);
    fwhm = param(3);
    normphots = sim_gaussian_photons(ngaussphots,mu,fwhm,ncurves);
    %% Simulate background photons using ITS
    nbckphots = param(4);
    bckphots = inverse_transform_sampling(tcal,bck,iOpt,nbckphots,ncurves);
    %% Simulate uniform photons
    noffsphots = param(5);
    offsphots = random("Uniform",0,T,noffsphots,ncurves);
    %% Put photons into big matrix
    photmatrix = vertcat(normphots,bckphots,offsphots);
    %% Bin photons into curves
    irf = bin_photons(tcal,nbins,ncurves,photmatrix);
    %% Apply Poisson weighting for photon counts
    irf(:) = random("Poisson",irf);
end

function phot = sim_gaussian_photons(nphots,mu,fwhm,ncurves)
    % Subroutine for simulating photons from Gaussian distribution

    sigma = fwhm / (2 * sqrt(2 * log(2)));
    phot = random("Normal",mu,sigma,nphots,ncurves);
end
