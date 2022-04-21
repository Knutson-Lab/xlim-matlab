function irf = shift_irf(irf,tcal,Q,iOpt)
    %SHIFT_IRF Subroutine for shifting IRF

    %% Housekeeping
    halfchan = tcal / 2;
    nbins = numel(irf);
    T = tcal * nbins;
    %% Prepare time axis to map IRF to
    tax = (halfchan:tcal:T)';
    itp = fit(tax,irf,iOpt);
    %% Scale Q (channels) to time
    tq = Q * tcal;
    irf = round(abs(itp(tax + tq)));
end

