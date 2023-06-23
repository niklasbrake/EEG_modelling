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
    xlim([0.1,40]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
subplot(1,3,2);
    plotwitherror(f,psd2,'CI','LineWidth',1,'color','r');
    hold on;
    plotwitherror(f,psd0,'CI','LineWidth',1,'color','k');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.1,40]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
subplot(1,3,3);
    plotwitherror(f,psd3,'CI','LineWidth',1,'color',[1,0,1]);
    hold on;
    plotwitherror(f,psd0,'CI','LineWidth',1,'color','k');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.1,40]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);


uni = load('E:\Research_Projects\004_Propofol\data\simulations\raw\peak_trend_sensitivity\trend_peak_interaction2.mat')

psd00 = mean(uni.P2(:,:,1),2);
psd11 = mean(uni.P2(:,:,6),2);
psd22 = mean(uni.P2(:,:,21),2);

figureNB(3,6.5);
subplot(3,1,1);
    plot(f,psd00,'color',[0.6,0.6,0.6],'LineWidth',1)
    hold on;
    plot(f,psd11,'color','k','LineWidth',1)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.1,100]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    xlabel('');
subplot(3,1,2);
    plot(f,psd00,'color',[0.6,0.6,0.6],'LineWidth',1)
    hold on;
    plot(f,psd22,'color','k','LineWidth',1)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.1,100]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    xlabel('');
subplot(3,1,3);
    plot(f,psd00,'color',[0.6,0.6,0.6],'LineWidth',1)
    hold on;
    plot(f,psd0,'LineWidth',1,'color','k');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.1,100]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);