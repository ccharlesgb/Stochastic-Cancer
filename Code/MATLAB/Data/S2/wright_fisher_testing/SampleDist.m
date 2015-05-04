function [] = SampleDist()
    figure();
    xAxis = 0:20
    
    %Normal
    subplot(121);
    data1 = load('wright_fisher_testing_PARAMS_K_21_N_1.0e9.0_d_100_u_1e-07_s_0.01.mat');
    SNames = fieldnames(data1); 
    disp(SNames)
    for i = 1:numel(SNames)
        if size(SNames{i}) < 6
            continue
        end
        if SNames{i}(1:6) == 'hist_n'
            semilogy(data1.hist_t, data1.(SNames{i}))
            hold all
        end
        %if SNames{i}
        %stuff = data1.(SNames{i});
        %disp(stuff)
    end
    xlabel('t')
    ylabel('n_j')
    subplot(122)
    for i = 1:numel(SNames)
        if size(SNames{i}) < 2
            continue
        end
        if (strcmp(SNames{i}(1:2),'t_') & strcmp(SNames{i}(end-2:end),'x_j'))
            semilogy(0:20, data1.(SNames{i}), 'o-')
            hold all
        end
    end
    xlabel('j')
    ylabel('n_j')
end

