function [] = fixtimefig()
    %Two Species Fix Time
    data_r1 = load('TwoSpeciesFix_Anal_sweepr1_N_10_SDP_5000.mat');
    data_u1 = load('TwoSpeciesFix_Anal_sweepu1_N_10_SDP_5000.mat');
    data_N = load('TwoSpeciesFix_Anal_sweepN_N_100_SDP_10000.mat');
    
    %Do r1 plot
    subplot(4,3,1);
    plot(data_r1.r1, data_r1.fixtime, '.k');
    hold on;
    plot(data_r1.r1, data_r1.num_fixtime);
    %xlabel('r1')
    ylabel('A. Fixation time')
    title('I')
    
    %Do u1 plot
    subplot(4,3,2);
    plot(data_u1.u1, data_u1.fixtime, '.k');
    hold on;
    plot(data_u1.u1, data_u1.num_fixtime);
    %xlabel('u1')
    %ylabel('Fixation time')
    title('II')
    
    %Do N plot
    subplot(4,3,3);
    plot(data_N.N, data_N.fixtime, '.k');
    hold on;
    plot(data_N.N, data_N.num_fixtime);
    %xlabel('N')
    %ylabel('Fixation time')
    title('III')
    
    %Deterministic Limit
    data_r1 = load('D:\Work\MPhys\Stochastic-Cancer\Code\MATLAB\Data\FixationTimeFig\ThreeSpeciesDeterm_sweepr1_N_10_SDP_5000.mat')
    data_u1 = load('D:\Work\MPhys\Stochastic-Cancer\Code\MATLAB\Data\FixationTimeFig\ThreeSpeciesDeterm_sweepu1_N_10_SDP_5000.mat')
    data_N = load('D:\Work\MPhys\Stochastic-Cancer\Code\MATLAB\Data\FixationTimeFig\ThreeSpeciesDeterm_sweepN_N_10_SDP_5000.mat')
    
    %Do r1 plot
    subplot(4,3,4);
    plot(data_r1.r1, data_r1.fixtime_r1, '.k');
    hold on;
    plot(data_r1.r1, data_r1.num_fixtime_r1);
    ylim([0,200])
    %xlabel('r1')
    ylabel('B. Fixation time')
    
    %Do u1 plot
    subplot(4,3,5);
    plot(data_u1.u1, data_u1.fixtime_u1, '.k');
    hold on;
    plot(data_u1.u1, data_u1.num_fixtime_u1);
    %xlabel('u1')
    %ylabel('Fixation time')
    
    %Do N plot
    subplot(4,3,6);
    plot(data_N.N, data_N.fixtime_N, '.k');
    hold on;
    plot(data_N.N, data_N.num_fixtime_N);
    %xlabel('N')
    %ylabel('Fixation time')
    
    
    %Naive three species
    data_r1 = load('ThreeSpeciesFix_Anal_sweepr1_N_10_SDP_5000.mat')
    data_u1 = load('ThreeSpeciesFix_Anal_sweepu1_N_10_SDP_5000.mat')
    data_N = load('ThreeSpeciesFix_Anal_sweepN_N_100_SDP_5000.mat')
    %Do r1 plot
    subplot(4,3,7);
    plot(data_r1.r1, data_r1.fixtime, '.k');
    hold on;
    plot(data_r1.r1, data_r1.num_fixtime);
    %xlabel('r1')
    ylabel('C. Fixation time')
    
    %Do u1 plot
    subplot(4,3,8);
    plot(data_u1.u1, data_u1.fixtime, '.k');
    hold on;
    plot(data_u1.u1, data_u1.num_fixtime);
    %xlabel('u1')
    %ylabel('Fixation time')
    
    %Do N plot
    subplot(4,3,9);
    plot(data_N.N, data_N.fixtime, '.k');
    hold on;
    plot(data_N.N, data_N.num_fixtime);
    %xlabel('N')
    %ylabel('Fixation time')
    
    %Systematic three species
    data_r1 = load('ThreeSpeciesSystematic_sweepr1_N_10_SDP_5000.mat')
    data_u1 = load('ThreeSpeciesSystematic_sweepu1_N_10_SDP_5000.mat')
    data_N =  load('ThreeSpeciesSystematic_sweepN_N_10_SDP_5000.mat')
    %Do r1 plot
    subplot(4,3,10);
    plot(data_r1.r1, data_r1.fixtime, '.k');
    hold on;
    plot(data_r1.r1, data_r1.num_fixtime);
    xlabel('r1')
    ylim([0,200])
    ylabel('D. Fixation time')
    
    %Do u1 plot
    subplot(4,3,11);
    plot(data_u1.u1, data_u1.fixtime, '.k');
    hold on;
    plot(data_u1.u1, data_u1.num_fixtime);
    xlabel('u1')
    %ylabel('Fixation time')
    
    %Do N plot
    subplot(4,3,12);
    plot(data_N.N, data_N.fixtime, '.k');
    hold on;
    plot(data_N.N, data_N.num_fixtime);
    xlabel('N')
    %ylabel('Fixation time')
    
end

