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


% uni = load('E:\Research_Projects\004_Propofol\data\simulations\raw\peak_trend_sensitivity\trend_peak_interaction2.mat')
uni = load('E:\Research_Projects\004_Propofol\data\simulations\raw\trend_peak_interaction\trend_peak_interaction');

psd00 = mean(uni.P2(:,:,1),2);
psd11 = mean(uni.P2(:,:,6),2);
psd22 = mean(uni.P2(:,:,21),2);
% psd000 = psd0./psd0(1000,:).*psd00(1000,:);

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
    ylim(10.^[-16.5,-14.5])
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
    ylim(10.^[-16.5,-14.5])
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
    ylim(10.^[-16.5,-14.5])