function [] = NewMethod()
    simCol = [0.0 0.0 1.0]
    oriCol = [0.0 0.5 0.0]
    tranCol = [1.0 0.0 0.0]
    negCol = [0.0 0.75 0.75]
    tran2Col = [1.0 0.0 0.0]
    
    minN = 1e6;
    maxN = 1e9;
    minU = 1e-8;
    maxU = 1e-5;
    minS = 1e-4;
    maxS = 1e-1;
    
    padding = 0.2;
    figure();
    dat_01 = load('WrightFisherSweepUNS_DPC_4_SDP_25_PARAM_K_21_N_1.0e9.0_d_100_u_1e-07_s_0.1_N_N.mat');
    %Sweep N
    subplot(1,3,1)
    semilogx(dat_01.NX, dat_01.Nt_20, 'o', 'Color', simCol);
    hold all
    semilogx(dat_01.NX, dat_01.Nt_20_TL, '^', 'Color', simCol);
    semilogx(dat_01.NX_orig, dat_01.Nt_20_orig, 'Color', oriCol);
    xlabel('N')
    ylabel('t_{20}')
    xlim([10.0 ^ (log10(minN) - padding), 10.0^(log10(maxN) + padding)]);
    
    %Sweep U
    subplot(1,3,2)
    loglog(dat_01.UX, dat_01.Ut_20, 'o', 'Color', simCol);
    hold all
    loglog(dat_01.UX_orig, dat_01.Ut_20_orig, 'Color', oriCol);
    xlabel('U')
    ylabel('t_{20}')
    xlim([10.0 ^ (log10(minU) - padding), 10.0^(log10(maxU) + padding)]);

    %Sweep S
    subplot(1,3,3)
    loglog(dat_01.SX, dat_01.St_20, 'o', 'Color', simCol);
    hold all
    loglog(dat_01.SX, dat_01.St_20_TL, '^', 'Color', simCol);
    xlim([10.0 ^ (log10(minS) - padding), 10.0^(log10(maxS) + padding)]);

    loglog(dat_01.SX_orig, dat_01.St_20_orig, 'Color', oriCol);
    xlabel('S');
    ylabel('t_{20}');
    legend('Simulation','Wright Fisher', 'Original');
    legend('boxoff');

end
