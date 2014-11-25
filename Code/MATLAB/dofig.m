function [] = dofig(name)
    set(gcf, 'Color', 'none');
    export_fig(strcat('MATLAB/',name), '-m1.5' ,'-transparent');
end

