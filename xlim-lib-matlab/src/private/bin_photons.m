function data = bin_photons(tcal,nbins,ncurves,phots)
    %BIN_PHOTONS Bins photons into discrete containers
    T = tcal * nbins;
    tbins = 0:tcal:T;
    type = class(phots);
    data = zeros(nbins,ncurves,type);
    data = fill_data(data,ncurves,phots,tbins);
end

function data = fill_data(data,ncurves,phots,tbins)
    for i = 1:ncurves
        %%% TODO: Any way to make this faster?
        data(:,i) = histcounts(phots(:,i),tbins);
    end
end
