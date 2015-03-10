function [] = DoTauHistory()
    figure();
    subplot(1,3,1)
    dat_01 = load('tauHistory_epsilon_0.1.mat');
    PlotIt(dat_01)
    
    subplot(1,3,2)
    dat_02 = load('tauHistory_epsilon_0.2.mat');
    PlotIt(dat_02)
    
    subplot(1,3,3)
    dat_03 = load('tauHistory_epsilon_0.3.mat');
    PlotIt(dat_03)
end

function [] = PlotIt(dat)
    hist(dat.tau)
    
    xlabel('Tau')
    ylabel('Frequency')
   % Get histogram patches
    ph = get(gca,'children');
    % Determine number of histogram patches
    N_patches = length(ph);
    for i = 1:N_patches
          % Get patch vertices
          vn = get(ph(i),'Vertices');
          % Adjust y location
          vn(:,2) = vn(:,2) + 1;
          % Reset data
          set(ph(i),'Vertices',vn)
    end
    % Change scale
    set(gca,'yscale','log')
end