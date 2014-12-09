function [] = import_cmap(x,y,z)
    z_flip = flipud(z);
    imagesc(x,y,z_flip);