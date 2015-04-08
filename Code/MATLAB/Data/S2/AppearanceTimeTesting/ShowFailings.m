function [] = ShowFailings()
    figure();
    xAxis = 0:20
    
    %Normal
    subplot(221);
    data1 = load('AppearanceTimeTesting_dpc_5_spd_10_params_K_21_N_1.0e9.0_d_100_u_1e-07_s_0.01.mat');
    plot(xAxis, data1.appear_sim, 'k.');
    hold all
    plot(xAxis, data1.appear_orig);
    title('(a) s = 1\times10^{-2}, u = 1\times10^{-7}');
    xlabel('Mutant type j');
    ylabel('Type j appearance time');
    
    %Low S
    subplot(222);
    data2 = load('AppearanceTimeTesting_dpc_5_spd_10_params_K_21_N_1.0e9.0_d_100_u_1e-07_s_0.0001.mat');
    plot(xAxis, data2.appear_sim, 'k.');
    hold all
    plot(xAxis, data2.appear_orig);
    title('(b) s = 1\times10^{-4}, u = 1\times10^{-7}');
    xlabel('Mutant type j');
    ylabel('Type j appearance time');
    
    %High U
    subplot(223);
    data3 = load('AppearanceTimeTesting_dpc_5_spd_10_params_K_21_N_1.0e9.0_d_100_u_1e-05_s_0.01.mat');
    plot(xAxis, data3.appear_sim, 'k.');
    hold all
    plot(xAxis, data3.appear_orig);
    title('(c) s = 1\times10^{-2}, u = 1\times10^{-5}');
    xlabel('Mutant type j');
    ylabel('Type j appearance time');
    
    %High U low S
    subplot(224);
    data4 = load('AppearanceTimeTesting_dpc_1_spd_10_params_K_21_N_1.0e9.0_d_100_u_1e-05_s_0.0001.mat');
    plot(xAxis, data4.appear_sim, 'k.');
    hold all
    plot(xAxis, data4.appear_orig);
    title('(d) s = 1\times10^{-4}, u = 1\times10^{-5}');
    xlabel('Mutant type j');
    ylabel('Type j appearance time');
end

