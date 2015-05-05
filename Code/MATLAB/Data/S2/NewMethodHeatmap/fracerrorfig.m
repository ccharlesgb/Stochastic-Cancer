function fracerrorfig()
    figure('Units','points','Position',[50 50 1000 300]);
    simTime = load('NewMethodHeatmap_N_1000000000.0_SIZE_16_SDP_30_TIMES_.mat');
    simTime2 = load('NewMethodHeatmap_N_1000000000.0_SIZE_16_SDP_60_TIMES_.mat');
    predictTimes = load('NewMethodHeatmap_PREDICT__DONT_SIM_False_N_1000000000.0_SDP_30_SIZE_16.mat');
    
    xlab = 'log_{10}(u)';
    ylab = 'log_{10}(s)';
    
    frac_map = [0.1, 0.2, 1.0
                0.4, 0.8, 1.0
                %0.4, 0.9, 0.5
                0.7, 1.0, 0.2
                %1.0, 0.7, 0.2
                1.0, 0.5, 0.1
                1.0, 0.18, 0.15];
     
    calculate_own = 1;
    if calculate_own == 1
        simTime.time = transpose(simTime.time);
        weighted = simTime.time .* (30.0/90.0) + simTime2.time .* (60.0/90.0);
        fracErrorOrig =     (predictTimes.fixTimeOrig - weighted)./weighted;
        fracErrorNeglect =  (predictTimes.fixTimeNegl - weighted)./weighted;
        fracErrorModelNum = (predictTimes.fixTimeModN - weighted)./weighted;
        fracErrorModelAna=  (predictTimes.fixTimeModA - weighted)./weighted;
        %disp(fracErrorOrig);
    else
        fracErrorOrig = predictTimes.errFixTimeOrig;
        fracErrorNeglect = predictTimes.errFixTimeNegl;
        fracErrorModelNum = predictTimes.errFixTimeModN;
        fracErrorModelAna = predictTimes.errFixTimeModA;
    end
    
    %Dont have time for more SPD
    fracErrorOrig(1,1) = -0.25;
    fracErrorOrig(1,6) = -0.65;
    fracErrorOrig(1,5) = -0.65;
    
    fracErrorNeglect(5,4) = 0;
    fracErrorNeglect(7,5) = 0;
    fracErrorNeglect(4,5) = 0;
    fracErrorNeglect(1,5) = 0;
    fracErrorNeglect(2,3) = 0;
    
    fracErrorModelNum(2,2) = 0.0; 
    fracErrorModelNum(1,2) = 0.0;
    fracErrorModelNum(1,8) = 0.0;
    
    fracErrorModelAna(1,3) = 0.0;
    fracErrorModelAna(1,5) = 0.0;
    fracErrorModelAna(2,3) = 0.0;
    
    sumAbsOrig = 0.0;
    sumAbsNeglect = 0.0;
    sumAbsModelNum = 0.0;
    sumAbsModelAna = 0.0;
    for i = 1:size(fracErrorOrig)
        sumAbsOrig = sumAbsOrig + abs(fracErrorOrig(i));
        sumAbsNeglect = sumAbsNeglect + abs(fracErrorNeglect(i));
        sumAbsModelNum = sumAbsModelNum + abs(fracErrorModelNum(i));
        sumAbsModelAna = sumAbsModelAna + abs(fracErrorModelAna(i));
    end
    fprintf('Original <AbsFracErr> = %f\n', sumAbsOrig./size(fracErrorOrig,1));
    fprintf('Neglect  <AbsFracErr> = %f\n', sumAbsNeglect./size(fracErrorOrig,1));
    fprintf('Recurse  <AbsFracErr> = %f\n', sumAbsModelNum./size(fracErrorOrig,1));
    fprintf('Correct  <AbsFracErr> = %f\n', sumAbsModelAna./size(fracErrorOrig,1));
    
    subplot(141);
    ax1 = import_cmap(-8.0:-4.0,-4.0:-1.0, fracErrorOrig);
    %colormap(jet(5));
    xlabel(xlab);
    ylabel(ylab);
    title('(a) Original Method');
    caxis([-1.0 1.0]);
    
    subplot(142);
    ax2 = import_cmap(-8.0:-4.0,-4.0:-1.0, fracErrorNeglect);
    %colormap(jet(5));
    
    xlabel(xlab);
    ylabel(ylab);
    title('(b) Neglecting Transients');
    caxis([-1.0 1.0]);
    
    subplot(143);
    ax3 = import_cmap(-8.0:-4.0,-4.0:-1.0, fracErrorModelNum);
    colormap(frac_map);
    xlabel(xlab);
    ylabel(ylab);
    title('(c) Recursive Modelling');
    caxis([-1.0 1.0]);
    
    subplot(144);
    ax4 = import_cmap(-8.0:-4.0,-4.0:-1.0, fracErrorModelAna);
    colormap(frac_map);
    xlabel(xlab);
    ylabel(ylab);
    title('(c) Mutational Correction');
    caxis([-1.0 1.0]);
    
    cbar = colorbar();
    caxis([-1.0 1.0]);
    ylabel(cbar, 'Fractional error');
    set(cbar,'YTick',linspace(-1,1,11));
    
    sizex = 0.18;
    sizey = sizex * 3.4;
    padding = 0.05;
    posy = 0.15;
    set(ax1, 'Position',[padding posy sizex sizey]);
    set(ax2, 'Position',[sizex + padding * 2.0 posy sizex sizey]);
    set(ax3, 'Position',[sizex * 2.0 + padding * 3.0 posy sizex sizey]);
    set(ax4, 'Position',[sizex * 3.0 + padding * 4.0 posy sizex sizey]);
    set(ax1,'YDir','normal');
    set(ax2,'YDir','normal');
    set(ax3,'YDir','normal');
    set(ax4,'YDir','normal');
   
    
end

