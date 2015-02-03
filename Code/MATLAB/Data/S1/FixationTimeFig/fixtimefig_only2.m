function [] = fixtimefig_only2()
    %Two Species Fix Time
    data_r1 = load('TwoSpeciesFix_Anal_sweepr1_N_10_SDP_5000.mat');
    data_u1 = load('TwoSpeciesFix_Anal_sweepu1_N_10_SDP_5000.mat');
    data_N = load('TwoSpeciesFix_Anal_sweepN_N_100_SDP_10000.mat');
    
    %Do r1 plot
    subplot(1,3,1);
    plot(data_r1.r1, data_r1.fixtime, '.k');
    hold on;
    plot(data_r1.r1, data_r1.num_fixtime);
    %xlabel('r1')
    ylabel('Fixation time')
    title('I')
    
    %Do u1 plot
    subplot(1,3,2);
    plot(data_u1.u1, data_u1.fixtime, '.k');
    hold on;
    plot(data_u1.u1, data_u1.num_fixtime);
    %xlabel('u1')
    %ylabel('Fixation time')
    title('II')
    
    %Do N plot
    subplot(1,3,3);
    plot(data_N.N, data_N.fixtime, '.k');
    hold on;
    plot(data_N.N, data_N.num_fixtime);
    %xlabel('N')
    %ylabel('Fixation time')
    title('III')
end

