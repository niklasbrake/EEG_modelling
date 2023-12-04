[full_model,synFun] = fittingmodel('lorenz');

load('E:\Research_Projects\004_Propofol\data\simulations\raw\mixed_input\trend_peak_interaction_mixed.mat')
f = 0.1:0.1:500;

psd0 = mean(P2(:,:,1),2);
psd1 = mean(P2(:,:,4),2);
psd2 = mean(P2(:,:,2),2);
psd3 = mean(P2(:,:,3),2);

figureNB(9,3)
subplot(1,3,1);
    plotwitherror(f,psd1,'CI','LineWidth',1,'color','b');
    hold on;
    plotwitherror(f,psd0,'CI','LineWidth',1,'color','k');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,50]);
    xticks([0.5,5,50]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    ylim(10.^[-16.5,-14])
subplot(1,3,2);
    plotwitherror(f,psd2,'CI','LineWidth',1,'color','r');
    hold on;
    plotwitherror(f,psd0,'CI','LineWidth',1,'color','k');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,50]);
    xticks([0.5,5,50]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    ylim(10.^[-16.5,-14])
subplot(1,3,3);
    plotwitherror(f,psd3,'CI','LineWidth',1,'color',[1,0,1]);
    hold on;
    plotwitherror(f,psd0,'CI','LineWidth',1,'color','k');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,50]);
    xticks([0.5,5,50]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    ylim(10.^[-16.5,-14])

gcaformat(gcf)


uni = load('E:\Research_Projects\004_Propofol\data\simulations\raw\trend_peak_interaction\trend_peak_interaction');

psd00 = mean(uni.P2(:,:,1),2);
psd11 = mean(uni.P2(:,:,6),2);
psd22 = mean(uni.P2(:,:,21),2);
f = 0.1:0.1:500;
% psd000 = psd0./psd0(1000,:).*psd00(1000,:);

figureNB(7.6,3);
axes('Position',[0.05, 0.23, 0.27, 0.70])
    plot(f,psd00,'color',[0.6,0.6,0.6],'LineWidth',1)
    hold on;
    plot(f,psd11,'color','k','LineWidth',1)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    ylim(10.^[-16.5,-14.5])
    xticks([0.5,5,50]);
axes('Position',[0.375, 0.23, 0.27, 0.70])
    plot(f,psd00,'color',[0.6,0.6,0.6],'LineWidth',1)
    hold on;
    plot(f,psd22,'color','k','LineWidth',1)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    ylim(10.^[-16.5,-14.5])
    xticks([0.5,5,50]);
axes('Position',[0.70, 0.23, 0.27, 0.70])
    plot(f,psd00,'color',[0.6,0.6,0.6],'LineWidth',1)
    hold on;
    plot(f,psd0,'LineWidth',1,'color','k');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    ylim(10.^[-16.5,-14.5])
    xticks([0.5,5,50]);
gcaformat(gcf)



%%%%%%% plot detrending figure
p0 = [0.13,0.0025,0.2,-0.85];
% psd0_detrend = 10*(log10(psd0./psd0(1))-synFun(f,p0)');
psd0_detrend = psd0./psd0(1)-10.^synFun(f,p0)';

p1 = [0.13,0.0025,0.25,-0.85];
% psd1_detrend = 10*(log10(psd1./psd1(1))-synFun(f,p1)');
psd1_detrend = psd1./psd1(1)-10.^synFun(f,p1)';

p2 = [0.16,0.0025,0.22,-1.27];
% psd2_detrend = 10*(log10(psd2./psd2(1))-synFun(f,p2)');
psd2_detrend = psd2./psd2(1)-10.^synFun(f,p2)';

p3 = [0.15,0.0025,0.25,-1.27];
% psd3_detrend = 10*(log10(psd3./psd3(1))-synFun(f,p3)');
psd3_detrend = psd3./psd3(1)-10.^synFun(f,p3)';


figureNB(14.5,6.5);
subplot(2,4,1);
    plot(f,psd0,'color','k','LineWidth',1);
    hold on;
    plot(f,psd0(1)*10.^synFun(f,p0),'color',[0.6,0.6,0.6],'LineStyle','-','LineWidth',1);
    % plot(f,psd0(1)*10.^synFun(f,p0),'color','k','LineStyle','--','LineWidth',1);
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
    plot(f,psd0,'color','k','LineWidth',1);
    hold on;
    plot(f,psd1,'color','b','LineWidth',1);
    hold on;
    plot(f,psd1(1)*10.^synFun(f,p1),'color',[0.6,0.6,0.6],'LineStyle','-','LineWidth',1);
    % plot(f,psd1(1)*10.^synFun(f,p1),'color','k','LineStyle','--','LineWidth',1);
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
    plot(f,psd0,'color','k','LineWidth',1);
    hold on;
    plot(f,psd2,'color','r','LineWidth',1);
    hold on;
    plot(f,psd2(1)*10.^synFun(f,p2),'color',[0.6,0.6,0.6],'LineStyle','-','LineWidth',1);
    % plot(f,psd2(1)*10.^synFun(f,p2),'color','k','LineStyle','--','LineWidth',1);
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
    plot(f,psd0,'color','k','LineWidth',1);
    hold on;
    plot(f,psd3,'color','m','LineWidth',1);
    hold on;
    plot(f,psd3(1)*10.^synFun(f,p3),'color',[0.6,0.6,0.6],'LineStyle','-','LineWidth',1);
    % plot(f,psd3(1)*10.^synFun(f,p3),'color','k','LineStyle','--','LineWidth',1);
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