function [] = dofig(name)
    linewidthnumber = 2;
    legendfontsize = 12;
    labelfontsize = 16;
    axisfontsize = 12;
    
    %custom_plot_handle_name = plot(x);
    %custom_legend_handle_name = legend();
    custom_x_label_name = get(gca, 'xlabel')
    custom_y_label_name = get(gca, 'ylabel')
    
    custom_title = get(gca, 'title');
    
    custom_legen_handle_name = legend(gca);
    
    
    %set (gcf, 'LineWidth', linewidthnumber);
    %set(custom_legen_handle_name,'FontSize',legendfontsize);
    %set(custom_x_label_name, 'FontSize', labelfontsize);
    %set(custom_y_label_name, 'FontSize', labelfontsize);
    
    %set(custom_title, 'FontSize', labelfontsize);
    
    %set(gca, 'FontSize', axisfontsize);
    set(gcf, 'Color', 'none');
    
    pause(0.1)
    
    export_fig(strcat('',name), '-m1.5' ,'-transparent', '-png', '-pdf', '-painters');
end

