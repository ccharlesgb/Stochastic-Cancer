function [] = fracerrorfig()
    figure()
    subplot(121)
    err_naive = load('NewMethodErrorMap_SDP=10_Size=8_N=1000000000.0.mat');

    import_cmap(0.5:1.5,0.5:1.5, err_naive.errFixTime1)

    colormap(jet(5))
    cbar = colorbar()
    caxis([-1.0 1.0])

    ylabel(cbar, 'Fractional error')
    
    subplot(122)

    import_cmap(0.5:1.5,0.5:1.5, err_naive.errFixTime2)

    colormap(jet(5))
    cbar = colorbar()
    caxis([-1.0 1.0])

    ylabel(cbar, 'Fractional error')
end

