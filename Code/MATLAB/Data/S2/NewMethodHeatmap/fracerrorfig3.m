function [] = fracerrorfig()
    figure('Units','points','Position',[50 50 800 300])
    simTime = load('NewMethodHeatmap_N_1000000000.0_SIZE_8_SDP_10_TIMES_.mat');
    predictTimes = load('NewMethodHeatmap_PREDICT__DONT_SIM_False_N_1000000000.0_SDP_10_SIZE_8.mat');
    
    xlab = 'log_{10}(s)';
    ylab = 'log_{10}(u)';
    
    frac_map = [0.1, 0.2, 1.0
                0.4, 0.8, 1.0
                0.5, 1.0, 0.3
                1.0, 0.5, 0.1
                1.0, 0.18, 0.15];
     
    calculate_own = 0 %Doesn't work oddly enough
    if calculate_own == 1
        fracErrorOrig =     (predictTimes.fixTimeOrig - simTime.time)./simTime.time;
        fracErrorNeglect =  (predictTimes.fixTimeNegl - simTime.time)./simTime.time;
        fracErrorModelNum = (predictTimes.fixTimeModN - simTime.time)./simTime.time;
        fracErrorModelAna=  (predictTimes.fixTimeModA - simTime.time)./simTime.time;
        disp(fracErrorOrig);
    else
        fracErrorOrig = predictTimes.errFixTimeOrig;
        fracErrorNeglect = predictTimes.errFixTimeNegl;
        fracErrorModelNum = predictTimes.errFixTimeModN;
        fracErrorModelAna = predictTimes.errFixTimeModA;
    end
    
    subplot(131)
    ax1 = import_cmap(-4.0:-1.0,-8.0:-5.0, fracErrorOrig);
    %colormap(jet(5));
    xlabel(xlab)
    ylabel(ylab)
    title('(a) Original Method');
    caxis([-1.0 1.0]);
    
    
    subplot(132);
    ax2 = import_cmap(-4.0:-1.0,-8.0:-5.0, fracErrorNeglect);
    %colormap(jet(5));
    
    xlabel(xlab)
    ylabel(ylab)
    title('(b) Neglecting Transients');
    caxis([-1.0 1.0]);
    
    subplot(133);
    ax3 = import_cmap(-4.0:-1.0,-8.0:-5.0, fracErrorModelAna);
    colormap(frac_map);
    xlabel(xlab);
    ylabel(ylab);
    title('(c) Modelling Transients');
    
    cbar = colorbar();
    caxis([-1.0 1.0]);
    ylabel(cbar, 'Fractional error');
    set(cbar,'YTick',linspace(-1,1,11));
    
    sizex = 0.24;
    sizey = sizex * 2.8;
    padding = 0.06;
    posy = 0.15;
    set(ax1, 'Position',[padding posy sizex sizey]);
    set(ax2, 'Position',[sizex + padding * 2.0 posy sizex sizey]);
    set(ax3, 'Position',[sizex * 2.0 + padding * 3.0 posy sizex sizey]);
    set(ax1,'YDir','normal');
    set(ax2,'YDir','normal');
    set(ax3,'YDir','normal');
   
    
end

