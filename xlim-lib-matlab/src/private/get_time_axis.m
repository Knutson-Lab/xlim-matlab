function tax = get_time_axis(tcal,trans)
    %GET_TIME_AXIS Returns time axis for curves
    [nbins] = size(trans);
    tst = tcal / 2;
    tend = tcal * nbins;
    ntax = numel(tcal);
    tax = zeros(nbins,ntax);
    tax = fill_tax(tax,tst,tcal,tend,ntax);
end

function tax = fill_tax(tax,tst,tcal,tend,ntax)
    for i = 1:ntax
        tax(:,i) = tst(i):tcal(i):tend(i);
    end
end
