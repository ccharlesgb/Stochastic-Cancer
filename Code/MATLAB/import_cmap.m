function [ax] = import_cmap(x,y,z)
    %z = flipud(z);
    ax = imagesc(x,y,z);
    ax = gca
    