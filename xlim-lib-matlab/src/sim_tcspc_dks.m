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
    % Check that number of irf/background curves match
    nirfbins = numel(irf);
    nbckbins = numel(bck);
    assert(nirfbins == nbckbins);
    % Check that parameter inputs align with specified number of
    % lifetimes/anisotropies
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

    % Check that all amplitude terms are in photon counts
    ntaucomp = 2 * ntau;
    ntaucounts = param(1:2:ntaucomp);
    nlincounts = vertcat(ntaucounts,param(end - 3:end));
    mustBeInteger(nlincounts);
    mustBeNonnegative(nlincounts);

    % Initialize photon matrix
    ntotalphots = sum(nlincounts);
    ctype = class(param);
    photonmatrix = zeros(ncurves,ntotalphots,ctype);
    bounds = cumsum(nlincounts);
    % Simulate photons from exponential distribution
    taus = param(2:2:ntaucomp);
    counter = 0;
    for i = 1:ntau
        localbound = bounds(i);
        localtau = taus(i);
        localtaucounts = ntaucounts(i);
        photonmatrix(:,(counter + 1):localbound) = ...
            random("Exponential",localtau,ncurves,localtaucounts);
        counter = counter + localbound;
    end
    % Simulate photons from IRF distribution using UTS
    halfchan = tcal / 2;
    nbins = numel(irf);
    T = tcal * nbins;
    tax = halfchan:tcal:T;
    itp = fit(tax,irf,iOpt);
    pdfirfitp = fit(tax,irf / integrate(itp,T,0),iOpt);
    % Shift IRF
    q = param(end - 3);
    pdfirfitp = fit(tax,pdfirfitp(tax + (tcal * q)));
    cdfitp = fit(tax,integrate(pdfirfitp,tax,0),iOpt);
    nscattercounts = param(end - 2);
    invcdfitp = fit(cdfitp(tax),tax,iOpt);
    localbound = bounds(end - 2);
    photonmatrix(:,(counter + 1):localbound) = ...
        invcdfitp(random("Uniform",0,1,ncurves,nscattercounts));
    counter = counter + localbound;
    % Repeat for background distribution using UTS
    itp = fit(tax,bck,iOpt);
    pdfitp = fit(tax,bck / integrate(itp,T,0),iOpt);
    cdfitp = fit(tax,integrate(pdfitp,tax,0),iOpt);
    nbackgroundcounts = param(end - 1);
    invcdfitp = fit(cdfitp(tax),tax,iOpt);
    localbound = bounds(end - 1);
    photonmatrix(:,(counter + 1):localbound) = ...
        invcdfitp(random("Uniform",0,1,ncurves,nbackgroundcounts));
    counter = counter + localbound;
    % Do for uniform distribution
    noffsetcounts = param(end);
    localbound = bounds(end);
    photonmatrix(:,(counter + 1):localbound) = ...
        invcdfitp(random("Uniform",0,T,ncurves,noffsetcounts));
    % Now package photons into bins
    tbins = 0:tcal:T;
    binnedPhots = zeros(nbins,ncurves,ctype);
    photonMatrixT = photonmatrix';
    for i = 1:ncurves
        binnedPhots(:,i) = histcounts(photonMatrixT(:,i),tbins);
    end
    % Convolve photons with PDF of IRF
    switch cMethod
        case "FFT"
            trans = round(abs(ifft(fftshift(fft(binnedPhots,[],1) .* fft(pdfirfitp(tax),[],1)),[],1)));
        otherwise
            error("Unsupported convolution!")
    end
end

