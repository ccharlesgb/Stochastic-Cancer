function [h, hcb] = imagescwithnan(a,cm,nanclr)
    % IMAGESC with NaNs assigning a specific color to NaNs
    a = flipud(a);
    %# find minimum and maximum
    amin=-0.5
    amax=0.5
    %# size of colormap
    n = size(cm,1);
    %# color step
    dmap=(amax-amin)/n;

    %# standard imagesc
    him = imagesc(0.5:1.5, 0.5:1.5, a);
    %# add nan color to colormap
    colormap([nanclr; cm]);
    %# changing color limits
    caxis([amin-dmap amax]);
    %# place a colorbar
    hcb = colorbar;
    %# change Y limit for colorbar to avoid showing NaN color
    ylim(hcb,[amin amax])
    ylabel(hcb, 'Fractional error')
    if nargout > 0
        h = him;
    end
    
        %# find minimum and maximum
    amin=0.0
    amax=0.5
    %# size of colormap
    n = size(jet(5),1);
    %# color step
    dmap=(amax-amin)/n;

    %# standard imagesc
    %# add nan color to colormap
    %olormap([nanclr; cm]);
    %# changing color limits
    caxis([amin-dmap amax]);