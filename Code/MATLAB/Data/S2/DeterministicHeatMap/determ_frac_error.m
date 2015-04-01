function [] = determ_frac_error()
    figure('Units','points','Position',[100 100 900 300])
    err_10 = load('DetermTauLeapFracError2_Size_20_SDP_60_N_10.0.mat');
    err_1000 = load('DetermTauLeapFracError2_Size_20_SDP_60_N_1000.0.mat');
    err_10000 = load('DetermTauLeapFracError2_Size_20_SDP_60_N_10000.0.mat');
    
    xlab = 'r_1';
    ylab = 'r_2';
    
    frac_map = [0.1, 0.2, 1.0
                0.4, 0.8, 1.0
                0.5, 1.0, 0.3
                1.0, 0.5, 0.1
                1.0, 0.18, 0.15];

    subplot(133);
    ax1 = import_cmap(0.5:1.5,0.5:1.5, err_10.fracError);
    colormap(frac_map);
    xlabel(xlab);
    ylabel(ylab);
    title('N = 10');
    caxis([-1.0 1.0]);
            
    subplot(131);
    ax2 = import_cmap(0.5:1.5,0.5:1.5, err_1000.fracError);
    %colormap(jet(5));
    xlabel(xlab)
    ylabel(ylab)
    title('N = 1000');
    caxis([-1.0 1.0]);
    
    subplot(132);
    ax3 = import_cmap(0.5:1.5,0.5:1.5, err_10000.fracError);
    %colormap(jet(5));
    
    xlabel(xlab)
    ylabel(ylab)
    title('N = 10000');
    caxis([-1.0 1.0]);
    
    colormap(frac_map);
    cbar = colorbar();
    
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

