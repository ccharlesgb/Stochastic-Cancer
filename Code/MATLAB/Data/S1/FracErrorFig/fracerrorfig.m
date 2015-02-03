function [] = fracerrorfig(type)
    if type == 1
        err_naive = load('D:\Work\MPhys\Stochastic-Cancer\Code\MATLAB\Data\FracErrorFig\three_species_error_Size_20_SDP_3000.mat');

        import_cmap(0.5:1.5,0.5:1.5, err_naive.frac_error)

        colormap(jet(5))
        cbar = colorbar()
        caxis([0.0,0.5])

        ylabel(cbar, 'Fractional error')
    elseif type == 2
        err_deter = load('D:\Work\MPhys\Stochastic-Cancer\Code\MATLAB\Data\FracErrorFig\three_species_errordeterm_Size_20_SDP_300.mat'); 
        imagescwithnan(err_deter.frac_error, jet(5), [0.5,0.5,0.5])
        
       % cbar = colorbar()
        
        %ylabel(cbar, 'Fractional Error')
    else
       print('invalid type') 
    end
end

