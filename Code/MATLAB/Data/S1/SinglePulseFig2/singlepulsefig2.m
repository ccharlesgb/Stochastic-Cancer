function [] = singlepulsefig
    figure();
    
    dat_1010 = load('DAT_PulseConstantArea_r1_1.0_r2_1.0_uquit_0.001_N_10_Dose_5.0.mat');
    param_1010 = load('PARAM_PulseConstantArea_r1_1.0_r2_1.0_uquit_0.001_N_10_Dose_5.0.mat');
    plotspdat(dat_1010, param_1010, 'Neutral' , 'r');

    dat_0910 = load('DAT_PulseConstantArea_r1_0.9_r2_1.0_uquit_0.001_N_10_Dose_5.0.mat');
    param_0910 = load('PARAM_PulseConstantArea_r1_0.9_r2_1.0_uquit_0.001_N_10_Dose_5.0.mat');
    plotspdat(dat_0910, param_0910, 'r1 = 0.9' ,'g');
    
    dat_1009 = load('DAT_PulseConstantArea_r1_1.0_r2_0.9_uquit_0.001_N_10_Dose_5.0.mat');
    param_1009 = load('PARAM_PulseConstantArea_r1_1.0_r2_0.9_uquit_0.001_N_10_Dose_5.0.mat');
    plotspdat(dat_1009, param_1009, 'r2 = 0.9' ,'b');
    
    dat_1110 = load('DAT_PulseConstantArea_r1_1.1_r2_1.0_uquit_0.001_N_10_Dose_5.0.mat');
    param_1110 = load('PARAM_PulseConstantArea_r1_1.1_r2_1.0_uquit_0.001_N_10_Dose_5.0.mat');
    plotspdat(dat_1110, param_1110, 'r1 = 1.1' ,'m');
    
    dat_1011 = load('DAT_PulseConstantArea_r1_1.0_r2_1.1_uquit_0.001_N_10_Dose_5.0.mat');
    param_1011 = load('PARAM_PulseConstantArea_r1_1.0_r2_1.1_uquit_0.001_N_10_Dose_5.0.mat');
    plotspdat(dat_1011, param_1011, 'r2 = 1.1' ,'c');

end

function plotspdat(fixes, params, label, col)
    subplot(2,1,1);
    
    %for i=1:5,
    %    PlotPulse(3.0, param_0910.u_pulse(i) / 0.001, 0.001, i);
    %    xlim([0,100]);
     %   ylabel('u(t)');
    %    xlabel('t');
    %end
    
    subplot(2,1,2);
    plot(fixes.PulseTime, fixes.FixTime - fixes.FixTime(1), col, 'DisplayName', label, 'Marker', '.');%,'.k');
    xlabel('t_{mut}');
    ylabel('Fixation Time');
    
    ratio = fixes.FixTime(6) / fixes.FixTime(1)
    
    hold on
end

function [] = PlotPulse(area, mut_factor, offset, colid)
    cols = ['r','g','b','m','c'];

    dataX = [];
    dataY = [];
    mutTime = area / (offset * (mut_factor - 1.0));
    dataX(1) = 0.0;
    dataX(2) = mutTime;
    dataX(3) = mutTime;
    dataX(4) = 150.0;
    
    dataY(1) = offset * mut_factor;
    dataY(2) = offset * mut_factor;
    dataY(3) = offset;
    dataY(4) = offset;

    plot(dataX, dataY, cols(colid));
    hold on
end

