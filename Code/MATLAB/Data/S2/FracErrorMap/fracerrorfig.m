function [] = fracerrorfig()
    figure('Units','points','Position',[100 100 900 300])
    err_naive = load('NewMethodErrorMap2_SDP=30_Size=8_N=1000000000.0.mat');
    err_naive.errFixTime1;
    
    xlab = 'log_{10}(s)';
    ylab = 'log_{10}(u)';
    
    frac_map = [0.1, 0.2, 1.0
                0.4, 0.8, 1.0
                0.5, 1.0, 0.3
                1.0, 0.5, 0.1
                1.0, 0.18, 0.15];

    
    subplot(131)
    ax1 = import_cmap(-4.0:-1.0,-8.0:-5.0, err_naive.errFixTime1)
    %colormap(jet(5));
    xlabel(xlab)
    ylabel(ylab)
    title('(a) Original Method');
    caxis([-1.0 1.0]);
    
    
    subplot(132);
    ax2 = import_cmap(-4.0:-1.0,-8.0:-5.0, err_naive.errFixTime2);
    %colormap(jet(5));
    
    xlabel(xlab)
    ylabel(ylab)
    title('(b) Neglecting Transients');
    caxis([-1.0 1.0]);
    
    subplot(133);
    ax3 = import_cmap(-4.0:-1.0,-8.0:-5.0, err_naive.errFixTime3);
    colormap(frac_map);
    xlabel(xlab);
    ylabel(ylab);
    title('(c) Modelling Transients');
    
    cbar = colorbar();
    caxis([-1.0 1.0]);
    ylabel(cbar, 'Fractional error');
    
    sizex = 0.25;
    sizey = sizex * 3.0;
    padding = 0.05;
    posy = 0.15;
    set(ax1, 'Position',[padding posy sizex sizey]);
    set(ax2, 'Position',[sizex + padding * 2.0 posy sizex sizey]);
    set(ax3, 'Position',[sizex * 2.0 + padding * 3.0 posy sizex sizey]);
    set(ax1,'YDir','normal');
    set(ax2,'YDir','normal');
    set(ax3,'YDir','normal');
   
    
    
end

