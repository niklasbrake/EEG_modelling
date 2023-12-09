load(fullfile(dataFolder,'simulations','trend_peak_interactions','parameter_changes'));
f = 0.1:0.1:500;

rhythms = load('E:\Research_Projects\004_Propofol\data\simulations\raw\peak_trend_sensitivity\rhythm_spectra.mat');
rhythms.psd = rhythms.psd([1,5,3,4,2]);


P3 = zeros(length(f),5,5);
P3(:,:,1) = mean(P2(:,:,1:5),2);
P3(:,:,2) = mean(P2(:,:,21:25),2);
P3(:,:,3) = mean(P2(:,:,11:15),2);
P3(:,:,4) = mean(P2(:,:,16:20),2);
P3(:,:,5) = mean(P2(:,:,6:10),2);


figureNB(4.5,4.8);
for i = 1:4
    subplot(2,2,i);
    plot(f,P3(:,1,5),'color',[0.6,0.6,1],'LineWidth',1);
    hold on;
    plot(f,P3(:,i+1,5),'color',[0,0,1],'LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.1,100])
    xticks([0.1,1,10,100]);
    xticklabels([0.1,1,10,100]);
    ylim([1e-17,1e-15]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    gcaformat;
end
figureNB(4.5,4.8);
for i = 1:4
    subplot(2,2,i);
    plot(f,P3(:,1,2),'color',[1,0.6,0.6],'LineWidth',1);
    hold on;
    plot(f,P3(:,i+1,2),'color',[1,0,0],'LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.1,100])
    xticks([0.1,1,10,100]);
    xticklabels([0.1,1,10,100]);
    ylim([1e-17,1e-15]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    gcaformat;
end



[full_model,synFun] = fittingmodel('eq1');
p(:,1) = [0.01,0.002,-0.34,-0.19];
p(:,2) = [0.03,0.002,-0.16,-0.41];
p(:,3) = [0.01,0.002,-0.08,-0.55];
p(:,4) = [0.01,0.0025,-1,-0.02];
p(:,5) = [0.011,0.0025,-0.25,-0.27];

m = 5;
k = 3;
ft1 = P3(1,1,1)*10.^synFun(f(:),p(:,1));
P1 = P3(:,1,m);
psd_detrend1 = 10*log10(P1./ft1);
ft2 = P3(1,k,1)*10.^synFun(f(:),p(:,k));
P2 = P3(:,k,m);
psd_detrend2 = 10*log10(P2./ft2);
pp0 = psd_detrend2-psd_detrend1;
pp1 = 10*log10(P2./P1);

R(1,2) = mean(pp0(and(f>=0.1,f<1)));
R(1,1) = mean(pp1(and(f>=0.1,f<1)));
R(2,2) = mean(pp0(and(f>=1,f<=3)));
R(2,1) = mean(pp1(and(f>=1,f<=3)));
R(3,2) = mean(pp0(and(f>=40,f<=60)));
R(3,1) = mean(pp1(and(f>=40,f<=60)));

x0 = 0.05;
figureNB(8.5,6.5);
axes('Position',[x0,0.6,0.24,0.33]);
    plot(f,P1,'color',[0.6,0.6,1],'LineWidth',1);
    hold on;
    plot(f,P2,'color','b','LineWidth',1);
    plot(f,ft1,'color','k','LineStyle','-','LineWidth',1);
    plot(f,ft2,'color','k','LineStyle','-','LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.5,50])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);
    ylim([10^(-16.5),10^(-15.5)])
    yticks([]);
    ylabel('log power');
      
axes('Position',[x0+1/3*1.03,0.6,0.24,0.33]);
    plot(f,psd_detrend2,'color','b','LineWidth',1);
    hold on;
    plot(f,psd_detrend1,'color',[0.6,0.6,1],'LineWidth',1);
    % ylim([-2,12]);
    ylabel('Detrended power (dB)');
    set(gca,'xscale','log')
    xlim([0.5,50])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);

axes('Position',[x0+2/3*1.03,0.6,0.24,0.33]);
    plot(f,pp0,'color','k','LineWidth',1);
    hold on;
    plot(f,pp1,'color',[0.6,0.6,0.6],'LineWidth',1);
    line([0.1,50],[0,0],'color','k','linestyle','--');
    % ylabel('\Delta detrended power (dB)')
    ylabel('\Delta power (dB)')
    set(gca,'xscale','log')
    xlim([0.5,50])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);
    ylim([-1,3]);

axes('Position',[x0,0.12,0.24,0.33]);
    plot(f,P1,'color',[0.6,0.6,1],'LineWidth',1);
    hold on;
    plot(f,P2,'color','b','LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.5,50])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);
    yticks([]);
    ylabel('log power');

axes('Position',[x0+1/3*1.03,0.12,0.24,0.33]);
    plot(f,pp1,'color','k','LineWidth',1);
    line([0.1,50],[0,0],'color','k','linestyle','--');
    ylabel('\Delta raw power (dB)')
    set(gca,'xscale','log')
    xlim([0.5,50])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);
    ylim([-1,3]);

axes('Position',[x0+2/3*1.03,0.12,0.24,0.33]);
    bar(R','EdgeColor','none')
    xticks([1,2]);
    ylim([-1.5,1.9]);
    xticklabels({'Raw','Detrended'})
    ylabel('\Delta power (dB)');

gcaformat(gcf)




[full_model,synFun] = fittingmodel('eq1');
p(:,1) = [0.01,0.002,-0.34,-0.19];
p(:,2) = [0.03,0.002,-0.16,-0.41];
p(:,3) = [0.01,0.002,-0.08,-0.55];
p(:,4) = [0.01,0.0025,-1,-0.02];
p(:,5) = [0.011,0.0025,-0.25,-0.27];

P3(4,:,2) = 0.75*P3(3,:,2)+0.25*P3(7,:,2);
P3(5,:,2) = 0.5*P3(3,:,2)+0.5*P3(7,:,2);
P3(6,:,2) = 0.25*P3(3,:,2)+0.75*P3(7,:,2);

m = 2;
k = 3;
ft1 = 2*P3(1,1,1)*10.^synFun(f(:),p(:,1));
P1 = P3(:,1,m)+P3(:,1,5);
psd_detrend1 = 10*log10(P1./ft1);
ft2 = 2*P3(1,k,1)*10.^synFun(f(:),p(:,k));
P2 = P3(:,k,m)+P3(:,k,5);
psd_detrend2 = 10*log10(P2./ft2);
pp0 = psd_detrend2-psd_detrend1;
pp1 = 10*log10(P2./P1);

R(1,2) = mean(pp0(and(f>=0.1,f<1)));
R(1,1) = mean(pp1(and(f>=0.1,f<1)));
R(2,2) = mean(pp0(and(f>=1,f<=3)));
R(2,1) = mean(pp1(and(f>=1,f<=3)));
R(3,2) = mean(pp0(and(f>=40,f<=60)));
R(3,1) = mean(pp1(and(f>=40,f<=60)));

x0 = 0.05;
figureNB(8.5,6.5);
axes('Position',[x0,0.6,0.24,0.33]);
    plot(f,P1,'color',[0.6,0.6,1],'LineWidth',1);
    hold on;
    plot(f,P2,'color','b','LineWidth',1);
    plot(f,ft1,'color','k','LineStyle','-','LineWidth',1);
    plot(f,ft2,'color','k','LineStyle','-','LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.1,50])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);
    ylim([10^(-16.5),10^(-15.5)]*2)
    yticks([]);
    ylabel('log power');

axes('Position',[x0+1/3*1.03,0.6,0.24,0.33]);
    plot(f,psd_detrend2,'color','b','LineWidth',1);
    hold on;
    plot(f,psd_detrend1,'color',[0.6,0.6,1],'LineWidth',1);
    % ylim([-2,12]);
    ylabel('Detrended power (dB)');
    set(gca,'xscale','log')
    xlim([0.1,50])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);

axes('Position',[x0+2/3*1.03,0.6,0.24,0.33]);
    plot(f,pp0,'color','k','LineWidth',1);
    hold on;
    plot(f,pp1,'color',[0.6,0.6,0.6],'LineWidth',1);
    line([0.1,50],[0,0],'color','k','linestyle','--');
    % ylabel('\Delta detrended power (dB)')
    ylabel('\Delta power (dB)')
    set(gca,'xscale','log')
    xlim([0.5,50])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);
    ylim([-1,3]);

axes('Position',[x0,0.12,0.24,0.33]);
    plot(f,P1,'color',[0.6,0.6,1],'LineWidth',1);
    hold on;
    plot(f,P2,'color','b','LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.5,50])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);
    yticks([]);
    ylabel('log power');

axes('Position',[x0+1/3*1.03,0.12,0.24,0.33]);
    plot(f,pp1,'color','k','LineWidth',1);
    line([0.1,50],[0,0],'color','k','linestyle','--');
    ylabel('\Delta raw power (dB)')
    set(gca,'xscale','log')
    xlim([0.5,50])
    xlabel('Frequency (Hz)');
    xticks([0.5,5,50]);
    ylim([-1,3]);

axes('Position',[x0+2/3*1.03,0.12,0.24,0.33]);
    bar(R','EdgeColor','none')
    xticks([1,2]);
    ylim([-1.5,1.9]);
    xticklabels({'Raw','Detrended'})
    ylabel('\Delta power (dB)');

gcaformat(gcf)