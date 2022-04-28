% function mpoint = get_mpoint(nparam,fxparam,gparam,ncurves)
%     %GET_MPOINT Returns pointer matrix for curve fitting
%     % Here, gparam needs to intersect with fxparam
%     gparam = gparam | fxparam;
%     Ml = false(nparam,ncurves); 
%     goffs = ncurves - 1;
%     if ncurves > 1
%         Ml(:,2:end) = repmat(gparam,1,goffs);
%     else
%         Ml(:) = gparam;
%     end
%     Mlx = ~Ml; mpoint = zeros(size(Ml));mpoint(Mlx) = 1:numel(find(Mlx));
%     if ncurves > 1
%         mpoint(Ml) = repmat(find(gparam),1,goffs);
%     else
%         mpoint(Ml) = find(gparam);
%     end
%     mpoint(fxparam,:) = -mpoint(fxparam,:);
% end

function mpoint = get_mpoint(nparam,fxparam,gparam,ncurves)
    % New approach with linear/matrix indexing
    gfxparam = gparam | fxparam;
    %% Generate coordinate map of curves
    [X,Y] = meshgrid(1:ncurves,1:nparam);
    %% Set any global and fixed parameters to the same index in X
    X(gfxparam,:) = repmat(X(gfxparam,1),1,ncurves);
    %% Create linear index map
    mpoint = sub2ind([nparam,ncurves],Y,X);
    %% Set any fixed values to 0
    mpoint(fxparam,:) = -mpoint(fxparam,:);
end
