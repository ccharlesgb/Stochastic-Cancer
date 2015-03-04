function [] = DoTauHistory()
    figure();
    subplot(1,3,1)
    dat_001 = load('tauHistory_epsilon_0.01.mat');
    PlotIt(dat_001)
    
    subplot(1,3,2)
    dat_01 = load('tauHistory_epsilon_0.1.mat');
    PlotIt(dat_01)
    
    subplot(1,3,3)
    dat_1 = load('tauHistory_epsilon_1.0.mat');
    PlotIt(dat_1)
end

function [] = PlotIt(dat)
    hist(dat.tau)
    hx = graph2d.constantline(dat.epsilon, 'LineStyle',':', 'Color',[.7 .7 .7]);
    changedependvar(hx,'x');
    
    xlabel('Tau')
    ylabel('Frequency')
end