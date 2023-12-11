[full_model,synFun] = fittingmodel('eq1');

load(fullfile(dataFolder,'simulations','trend_peak_interactions','mixed_dynamics'));

psd_baseline = mean(psd_baseline,2);
psd_high_oscillation = mean(psd_high_oscillation,2);
psd_high_aperiodic = mean(psd_high_aperiodic,2);
psd_high_both = mean(psd_high_both,2);
psd_asynch = mean(psd_asynch,2);
psd_aperiodic_only = mean(psd_aperiodic_only,2);
psd_oscillation_only = mean(psd_oscillation_only,2);

figureNB(7.6,3);
axes('Position',[0.05, 0.23, 0.27, 0.70])
    plot(f,psd_asynch,'color',[0.6,0.6,0.6],'LineWidth',1)
    hold on;
    plot(f,psd_aperiodic_only,'color','k','LineWidth',1)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    ylim(10.^[-16.5,-14])
    xticks([0.5,5,50]);
axes('Position',[0.375, 0.23, 0.27, 0.70])
    plot(f,psd_asynch,'color',[0.6,0.6,0.6],'LineWidth',1)
    hold on;
    plot(f,psd_oscillation_only,'color','k','LineWidth',1)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    ylim(10.^[-16.5,-14])
    xticks([0.5,5,50]);
axes('Position',[0.70, 0.23, 0.27, 0.70])
    plot(f,psd_asynch,'color',[0.6,0.6,0.6],'LineWidth',1)
    hold on;
    plot(f,psd_baseline,'LineWidth',1,'color','k');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    ylim(10.^[-16.5,-14])
    xticks([0.5,5,50]);
gcaformat(gcf)



%%%%%%% plot detrending figure
p0 = [0.13,0.0025,0.2,-0.85];
psd0_detrend = 10*(log10(psd_baseline./psd_baseline(1))-synFun(f,p0)');
% psd0_detrend = psd_baseline./psd_baseline(1)-10.^synFun(f,p0)';

p1 = [0.13,0.0025,0.25,-0.85];
psd1_detrend = 10*(log10(psd_high_oscillation./psd_high_oscillation(1))-synFun(f,p1)');
% psd1_detrend = psd_high_oscillation./psd_high_oscillation(1)-10.^synFun(f,p1)';

p2 = [0.16,0.0025,0.22,-1.27];
psd2_detrend = 10*(log10(psd_high_aperiodic./psd_high_aperiodic(1))-synFun(f,p2)');
% psd2_detrend = psd_high_aperiodic./psd_high_aperiodic(1)-10.^synFun(f,p2)';

p3 = [0.15,0.0025,0.25,-1.27];
psd3_detrend = 10*(log10(psd_high_both./psd_high_both(1))-synFun(f,p3)');
% psd3_detrend = psd_high_both./psd_high_both(1)-10.^synFun(f,p3)';


figureNB(14.5,6.5);
subplot(2,4,1);
    plot(f,psd_baseline,'color','k','LineWidth',1);
    hold on;
    plot(f,psd_baseline(1)*10.^synFun(f,p0),'color',[0.6,0.6,0.6],'LineStyle','-','LineWidth',1);
    % plot(f,psd_baseline(1)*10.^synFun(f,p0),'color','k','LineStyle','--','LineWidth',1);
    set(gca,'yscale','log');
    ylim(10.^[-16.5,-14])
    yticks([]);
    ylabel('log power');
    set(gca,'xscale','log');
    xlim([0.5,100])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);
subplot(2,4,5);
    plot(f,psd0_detrend,'color','k','LineWidth',1);
    ylim([-2,12]);
    ylabel('Detrended power (dB)');
    set(gca,'xscale','log');
    xlim([0.5,100])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);

subplot(2,4,2);
    plot(f,psd_baseline,'color','k','LineWidth',1);
    hold on;
    plot(f,psd_high_oscillation,'color','b','LineWidth',1);
    hold on;
    plot(f,psd_high_oscillation(1)*10.^synFun(f,p1),'color',[0.6,0.6,0.6],'LineStyle','-','LineWidth',1);
    % plot(f,psd_high_oscillation(1)*10.^synFun(f,p1),'color','k','LineStyle','--','LineWidth',1);
    set(gca,'yscale','log');
    ylim(10.^[-16.5,-14])
    yticks([]);
    ylabel('log power');
    set(gca,'xscale','log');
    xlim([0.5,100])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);
subplot(2,4,6);
    plot(f,psd0_detrend,'color','k','LineWidth',1);
    ylim([-2,12]);
    ylabel('Detrended power (dB)');
    hold on;
    plot(f,psd1_detrend,'color','b','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.5,100])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);

subplot(2,4,3);
    plot(f,psd_baseline,'color','k','LineWidth',1);
    hold on;
    plot(f,psd_high_aperiodic,'color','r','LineWidth',1);
    hold on;
    plot(f,psd_high_aperiodic(1)*10.^synFun(f,p2),'color',[0.6,0.6,0.6],'LineStyle','-','LineWidth',1);
    % plot(f,psd_high_aperiodic(1)*10.^synFun(f,p2),'color','k','LineStyle','--','LineWidth',1);
    set(gca,'yscale','log');
    ylim(10.^[-16.5,-14])
    yticks([]);
    ylabel('log power');
    set(gca,'xscale','log');
    xlim([0.5,100])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);
subplot(2,4,7);
    plot(f,psd0_detrend,'color','k','LineWidth',1);
    ylim([-2,12]);
    ylabel('Detrended power (dB)');
    hold on;
    plot(f,psd2_detrend,'color','r','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.5,100])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);


subplot(2,4,4);
    plot(f,psd_baseline,'color','k','LineWidth',1);
    hold on;
    plot(f,psd_high_both,'color','m','LineWidth',1);
    hold on;
    plot(f,psd_high_both(1)*10.^synFun(f,p3),'color',[0.6,0.6,0.6],'LineStyle','-','LineWidth',1);
    % plot(f,psd_high_both(1)*10.^synFun(f,p3),'color','k','LineStyle','--','LineWidth',1);
    set(gca,'yscale','log');
    ylim(10.^[-16.5,-14])
    yticks([]);
    ylabel('log power');
    set(gca,'xscale','log');
    xlim([0.5,100])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);
subplot(2,4,8);
    plot(f,psd0_detrend,'color','k','LineWidth',1);
    ylim([-2,12]);
    ylabel('Detrended power (dB)');
    hold on;
    plot(f,psd3_detrend,'color','m','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.5,100])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);
gcaformat(gcf)