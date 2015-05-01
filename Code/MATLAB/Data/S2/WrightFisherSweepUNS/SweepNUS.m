function [] = NewMethod()
    simCol = [0.0 0.0 1.0]
    oriCol = [0.0 0.5 0.0]
    tranCol = [1.0 0.0 0.0]
    negCol = [0.0 0.75 0.75]
    tran2Col = [1.0 0.0 0.0]

    figure();
    dat_01 = load('WrightFisherSweepUNS_DPC_4_SDP_25_PARAM_K_21_N_1.0e9.0_d_100_u_1e-07_s_0.1_N_N.mat');
    %Sweep N
    subplot(1,3,1)
    loglog(dat_01.NX, dat_01.Nt_20, 'o', 'Color', simCol);
    hold all
    loglog(dat_01.NX, dat_01.Nt_20_TL, '^', 'Color', simCol);
    loglog(dat_01.NX_orig, dat_01.Nt_20_orig, 'Color', oriCol);
    xlabel('N')
    ylabel('t_{20}')
    
    %Sweep U
    subplot(1,3,2)
    loglog(dat_01.UX, dat_01.Ut_20, 'o', 'Color', simCol);
    hold all
    loglog(dat_01.UX_orig, dat_01.Ut_20_orig, 'Color', oriCol);
    xlabel('U')
    ylabel('t_{20}')

    %Sweep S
    subplot(1,3,3)
    loglog(dat_01.SX, dat_01.St_20, 'o', 'Color', simCol);
    hold all
    loglog(dat_01.SX, dat_01.St_20_TL, '^', 'Color', simCol);

    loglog(dat_01.SX_orig, dat_01.St_20_orig, 'Color', oriCol);
    xlabel('S');
    ylabel('t_{20}');
    legend('Simulation','Wright Fisher', 'Original');
    legend('boxoff');

end
