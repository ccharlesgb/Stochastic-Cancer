function [] = VarUPlots()
    figure();

    %SIN 3
    subplot(231);
    data1 = load('VARU_APPTIME_SIN3.mat');
    plt_title = '(a) u=10^{-7}, 10^{-10}, 10^{-7}, 10^{-4}...';
    PlotVarU(data1, plt_title);
    %SIN 4
    subplot(234);
    data1 = load('VARU_APPTIME_SIN4.mat');
    plt_title = '(b) u=10^{-7}, 10^{-11}, 10^{-7}, 10^{-3}...';
    PlotVarU(data1, plt_title);
    
    %EXP INC
    subplot(232);
    data1 = load('VARU_APPTIME_EXP_INC.mat');
    plt_title = '(c) Exponential increase 10^{-9} to 10^{-6}';
    PlotVarU(data1, plt_title);
    %EXP DEC
    subplot(235);
    data1 = load('VARU_APPTIME_EXP_DEC.mat');
    plt_title = '(d) Exponential decrease 10^{-9} to 10^{-6}';
    PlotVarU(data1, plt_title);
    
    %SWITCH INC
    subplot(233);
    data1 = load('VARU_APPTIME_SWI_INC.mat');
    plt_title = '(e) Instant increase 10^{-9} to 10^{-6}';
    PlotVarU(data1, plt_title);
    %SWITCH DEC
    subplot(236);
    data1 = load('VARU_APPTIME_SWI_DEC.mat');
    plt_title = '(f) Instant decrease 10^{-9} to 10^{-6}';
    PlotVarU(data1, plt_title);
    
    legend('Simulation', 'Neglect', 'Correct', 'Recursive','Location','northwest')
    legend('boxoff');
end

function [] = PlotVarU(data, plt_title)
    xAxis = 0:20
    plot(xAxis, data.appear_sim, 'k.');
    hold all
    plot(xAxis, data.appear_neglect);
    plot(xAxis, data.appear_correct);
    plot(xAxis, data.appear_recursive);
    title(plt_title);
    xlabel('Mutant type j');
    ylabel('Type j appearance time');
end
