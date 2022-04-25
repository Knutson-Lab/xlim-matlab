function mpoint = get_mpoint(nparam,fxparam,gparam,ncurves)
    %GET_MPOINT Returns pointer matrix for curve fitting
    % Here, gparam needs to intersect with fxparam
    gparam = gparam | fxparam;
    Ml = false(nparam,ncurves); 
    goffs = ncurves - 1;
    if ncurves > 1
        Ml(:,2:end) = repmat(gparam,1,goffs);
    else
        Ml(:) = gparam;
    end
    Mlx = ~Ml; mpoint = zeros(size(Ml));mpoint(Mlx) = 1:numel(find(Mlx));
    mpoint(Ml) = repmat(find(gparam),1,goffs);
    mpoint(fxparam,:) = -mpoint(fxparam,:);
end

