function phots = inverse_transform_sampling(tcal,curve,iOpt,nphot,ncurves)
    %INVERSE_TRANSFORM_SAMPLING Subroutine that performs inverse transform
    %sampling from custom distribution
    ncounts = sum(curve);
    if ncounts > 0
        % Prepare time axis to plot with curve
        tcal2 = tcal / 2;
        nbins = numel(curve);
        T = nbins * tcal;
        tax = (tcal2:tcal:T)';
        % Generate inverse CDF
        itp = fit(tax,curve,iOpt);
        pdfitp = fit(tax,curve / integrate(itp,T,0),iOpt);
        cdfitp = fit(tax,integrate(pdfitp,tax,0),iOpt);
        invcdfitp = fit(cdfitp(tax),tax,iOpt);
        % Assign uniform photons to inverted distribution function
        phots = reshape(invcdfitp(random("Uniform",0,1,nphot,ncurves)),...
            nphot,ncurves);
    else
        % Return empty matrix
        phots = [];
    end
end

