function [] = NewMethod()
    figure();
    %Sweep N
    subplot(2,3,1)
    dat_01 = load('SweepN.mat');
    semilogx(dat_01.N, dat_01.Nt_20, 'o');
    hold all
    semilogx(dat_01.N, dat_01.Nt_20_anal1, '--');
    semilogx(dat_01.N, dat_01.Nt_20_anal2_transient);
    xlabel('N')
    ylabel('t_{20}')
    %Error N
    subplot(2,3,4)
    semilogx(dat_01.N, (dat_01.Nt_20-dat_01.Nt_20_anal1)./dat_01.Nt_20, '--g');
    hold all;
    semilogx(dat_01.N, (dat_01.Nt_20-dat_01.Nt_20_anal2_transient)./dat_01.Nt_20, 'r');
    hy = graph2d.constantline(0, 'Color',[.7 .7 .7]);
    changedependvar(hy,'y');
    xlim([dat_01.N(1) dat_01.N(length(dat_01.N))])
    xlabel('N')
    ylabel('Error')
    
    %Sweep U
    subplot(2,3,2)
    dat_02 = load('SweepU.mat');
    semilogx(dat_02.U, dat_02.Ut_20, 'o');
    hold all
    semilogx(dat_02.U, dat_02.Ut_20_anal1, '--');
    semilogx(dat_02.U, dat_02.Ut_20_anal2_transient);
    xlabel('U')
    ylabel('t_{20}')
    %Error U
    subplot(2,3,5)
    semilogx(dat_02.U, (dat_02.Ut_20-dat_02.Ut_20_anal1)./dat_02.Ut_20, '--g');
    hold all;
    semilogx(dat_02.U, (dat_02.Ut_20-dat_02.Ut_20_anal2_transient)./dat_02.Ut_20, 'r');
    hy = graph2d.constantline(0, 'Color',[.7 .7 .7]);
    changedependvar(hy,'y');
    xlim([dat_02.U(1) dat_02.U(length(dat_02.U))])
    xlabel('U')
    ylabel('Error')

    %Sweep S
    subplot(2,3,3)
    dat_03 = load('SweepS.mat');
    semilogx(dat_03.S, dat_03.St_20 ,'o');
    hold all
    semilogx(dat_03.S, dat_03.St_20_anal1, '--');
    semilogx(dat_03.S, dat_03.St_20_anal2_transient);
    xlabel('S')
    ylabel('t_{20}')
    %Error S
    subplot(2,3,6)
    semilogx(dat_03.S, (dat_03.St_20-dat_03.St_20_anal1)./(dat_03.St_20), '--g');
    hold all;
    semilogx(dat_03.S, (dat_03.St_20-dat_03.St_20_anal2_transient)./(dat_03.St_20), 'r');
    hy = graph2d.constantline(0, 'Color',[.7 .7 .7]);
    changedependvar(hy,'y');
    xlim([dat_03.S(1) dat_03.S(length(dat_03.S))])
    xlabel('S')
    ylabel('Error')
end
