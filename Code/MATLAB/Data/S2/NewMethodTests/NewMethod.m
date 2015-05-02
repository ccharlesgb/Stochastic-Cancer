function [] = NewMethod()
    simCol = [0.0 0.0 1.0]
    oriCol = [0.0 0.5 0.0]
    tranCol = [1.0 0.0 0.0]
    negCol = [0.0 0.75 0.75]
    tran2Col = [0.75 0.25 0.75]

    figure();
    %Sweep N
    subplot(2,3,1)
    dat_01 = load('SweepN.mat');
    semilogx(dat_01.N, dat_01.Nt_20, 'o', 'Color', simCol);
    hold all
    semilogx(dat_01.N, dat_01.Nt_20_anal1, 'Color', oriCol);
    semilogx(dat_01.N, dat_01.Nt_20_anal2_recursive, 'Color', tranCol);
    semilogx(dat_01.N, dat_01.Nt_20_anal4_correct, 'Color', tran2Col);
    semilogx(dat_01.N, dat_01.Nt_20_anal3_neglect, 'Color', negCol);
    xlabel('N')
    ylabel('t_{20}')
    %Error N
    subplot(2,3,4)
    semilogx(dat_01.N, (dat_01.Nt_20-dat_01.Nt_20_anal1)./dat_01.Nt_20, 'Color', oriCol);
    hold all;
    semilogx(dat_01.N, abs(dat_01.Nt_20-dat_01.Nt_20_anal2_recursive)./dat_01.Nt_20, 'Color', tranCol);
    semilogx(dat_01.N, abs(dat_01.Nt_20-dat_01.Nt_20_anal4_correct)./dat_01.Nt_20, 'Color', tran2Col);
    semilogx(dat_01.N, abs(dat_01.Nt_20-dat_01.Nt_20_anal3_neglect)./dat_01.Nt_20, 'Color', negCol);
    hy = graph2d.constantline(0, 'Color',[.7 .7 .7]);
    changedependvar(hy,'y');
    xlim([dat_01.N(1) dat_01.N(length(dat_01.N))])
    xlabel('N')
    ylabel('Absolute Fractional Error')
    ylim([0 1.0])
    
    %Sweep U
    subplot(2,3,2)
    dat_02 = load('SweepU.mat');
    semilogx(dat_02.U, dat_02.Ut_20, 'o', 'Color', simCol);
    hold all
    semilogx(dat_02.U, dat_02.Ut_20_anal1, 'Color', oriCol);
    semilogx(dat_02.U, dat_02.Ut_20_anal2_recursive, 'Color', tranCol);
    semilogx(dat_02.U, dat_02.Ut_20_anal4_correct, 'Color', tran2Col);
    semilogx(dat_02.U, dat_02.Ut_20_anal3_neglect, 'Color', negCol);
    xlabel('U')
    ylabel('t_{20}')
    %Error U
    subplot(2,3,5)
    semilogx(dat_02.U, (dat_02.Ut_20-dat_02.Ut_20_anal1)./dat_02.Ut_20, 'Color', oriCol);
    hold all;
    semilogx(dat_02.U, abs(dat_02.Ut_20-dat_02.Ut_20_anal2_recursive)./dat_02.Ut_20, 'Color', tranCol);
    semilogx(dat_02.U, abs(dat_02.Ut_20-dat_02.Ut_20_anal4_correct)./dat_02.Ut_20, 'Color', tran2Col);
    semilogx(dat_02.U, abs(dat_02.Ut_20-dat_02.Ut_20_anal3_neglect)./dat_02.Ut_20, 'Color', negCol);
    hy = graph2d.constantline(0, 'Color',[.7 .7 .7]);
    changedependvar(hy,'y');
    xlim([dat_02.U(1) dat_02.U(length(dat_02.U))])
    xlabel('U')
    ylabel('Absolute Fractional Error')
    ylim([0 1.0])

    %Sweep S
    subplot(2,3,3)
    dat_03 = load('SweepS.mat');
    semilogx(dat_03.S, dat_03.St_20 ,'o', 'Color', simCol);
    hold all
    semilogx(dat_03.S, dat_03.St_20_anal1, 'Color', oriCol);
    semilogx(dat_03.S, dat_03.St_20_anal2_recursive, 'Color', tranCol);
    semilogx(dat_03.S, dat_03.St_20_anal4_correct, 'Color', tran2Col);
    semilogx(dat_03.S, dat_03.St_20_anal3_neglect, 'Color', negCol);
    xlabel('S')
    ylabel('t_{20}')
    legend('Simulation', 'Original', 'Recursive', 'Single Correction', 'Neglect');
    legend('boxoff');
    %Error S
    ax = subplot(2,3,6);
    semilogx(dat_03.S, (dat_03.St_20-dat_03.St_20_anal1)./(dat_03.St_20), 'Color', oriCol);
    hold all;
    semilogx(dat_03.S, abs(dat_03.St_20-dat_03.St_20_anal2_recursive)./(dat_03.St_20), 'Color', tranCol);
    semilogx(dat_03.S, abs(dat_03.St_20-dat_03.St_20_anal4_correct)./(dat_03.St_20), 'Color', tran2Col);
    semilogx(dat_03.S, abs(dat_03.St_20-dat_03.St_20_anal3_neglect)./(dat_03.St_20), 'Color', negCol);
    hy = graph2d.constantline(0, 'Color',[.7 .7 .7]);
    changedependvar(hy,'y');
    xlim([dat_03.S(1) dat_03.S(length(dat_03.S))]);
    xlabel('S');
    ylabel('Absolute Fractional Error');
    ylim([0 1.0])
    
    legend('Original', 'Recursive', 'Single Correction', 'Neglect');
    legend('boxoff')
end
