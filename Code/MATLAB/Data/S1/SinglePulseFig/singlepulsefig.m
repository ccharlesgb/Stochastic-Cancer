function [] = singlepulsefig

    dat_3 = load('DAT_PulseConstantArea_r1_1.0_r2_1.0_uquit_0.001_N_10_Dose_3.0.mat');
    param_3 = load('PARAM_PulseConstantArea_r1_1.0_r2_1.0_uquit_0.001_N_10_Dose_3.0.mat')
    dat_5 = load('DAT_PulseConstantArea_r1_1.0_r2_1.0_uquit_0.001_N_10_Dose_5.0.mat');
    param_5 = load('PARAM_PulseConstantArea_r1_1.0_r2_1.0_uquit_0.001_N_10_Dose_5.0.mat')
    
    subplot(2,1,1);
  
    for i=1:5,
        PlotPulse(3.0, param_3.u_pulse(i) / 0.001, 0.001, i);
        xlim([0,100]);
        ylabel('u(t)');
        xlabel('t');
    end
    
    subplot(2,1,2);
    plot(dat_3.PulseTime, dat_3.FixTime)%,'.k');
    xlabel('t_{mut}');
    ylabel('Fixation Time');
    
    figure();
    subplot(2,1,1);
  
    for i=1:5,
        PlotPulse(5.0, param_5.u_pulse(i) / 0.001, 0.001, i);
        xlim([0,100]);
        ylabel('u(t)');
        xlabel('t');
    end

    subplot(2,1,2);
    plot(dat_5.PulseTime, dat_5.FixTime)%,'.k');
    xlabel('t_{mut}');
    ylabel('Fixation Time');
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

