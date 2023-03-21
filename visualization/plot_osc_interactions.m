red = clrsPT.qualitative_CM.red*0.85;
blue = clrsPT.qualitative_CM.blue*0.85;
grey = clrsPT.qualitative_CM.grey*0.75;
black = [0,0,0];


load('C:\Users\brake\Documents\temp\analzye_simulations_10.mat')
% P1(:,:,7:10) = nan;
P_tau_crit = P(:,:,2);
P_crit = P(:,:,4);
P_baseline = 1.2*P(:,:,6);
P_baseline_osc = P(:,:,7);
P_tau_osc = P(:,:,8);
P_tau = P(:,:,9);
P_crit_osc = P(:,:,10);

figureNB(11.5,8);
ax = axes('Position',[0.08,0.65,0.19,0.3]);
labelpanel(ax.Position(1)-0.05,ax.Position(2)+ax.Position(4)-0.01,'b',true);
    plot(f,mean(P_baseline_osc,2),'LineWidth',1,'color','k')
    hold on;
    plot(f,mean(P_crit_osc,2),'LineWidth',1,'color',red)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.1,100])
    ylim(10.^[-17,-14]);
    xticks(10.^[-1:2]);
    xticklabels([0.1,1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel('log power')
    % ylabel(['PSD (' char(956) 'V^2/Hz)']);
    text(0.13,0.8*4e-17,'m=0','color','k','FontSize',6)
    text(0.13,0.8*2e-17,'m=0.98','color',red,'FontSize',6)
    yticks([]);

    title('Single network','FontWeight','normal','FontSize',7);
ax = axes('Position',[0.42,0.65,0.19,0.3]);
labelpanel(ax.Position(1)-0.05,ax.Position(2)+ax.Position(4)-0.01,'c',true);
    plot(f,mean(P_baseline_osc+P_baseline,2),'LineWidth',1,'color','k')
    hold on;
    plot(f,mean(P_baseline_osc+P_crit,2),'LineWidth',1,'color',red)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.1,100])
    ylim(10.^[-17,-14]);
    xticks(10.^[-1:2]);
    xticklabels([0.1,1,10,100]);
    xlabel('Frequency (Hz)');
    text(0.13,0.8*4e-17,'m_1=0 / m_2=0','color','k','FontSize',6)
    text(0.13,0.8*2e-17,'m_1=0 / m_2=0.98','color',red,'FontSize',6)
    ylabel('log power')
    % ylabel(['PSD (' char(956) 'V^2/Hz)']);
    yticks([]);
    title('Mixed network','FontWeight','normal','FontSize',7);
ax = axes('Position',[0.76,0.65,0.19,0.3]);
labelpanel(ax.Position(1)-0.05,ax.Position(2)+ax.Position(4)-0.01,'d',true);
    plot(f,mean(P_crit+P_baseline_osc,2),'LineWidth',1,'color',red)
    hold on;
    plot(f,mean(P_tau_crit+P_tau_osc,2),'LineWidth',1,'color',blue)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.1,100])
    ylim(10.^[-16.5,-13.5]);
    xticks(10.^[-1:2]);
    xticklabels([0.1,1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel('log power')
    % ylabel(['PSD (' char(956) 'V^2/Hz)']);
    yticks([]);
    text(0.13,0.8*4*10^-16.5,'\tau =10 ms','color',red,'FontSize',6)
    text(0.13,0.8*2*10^-16.5,'\tau =30 ms','color',blue,'FontSize',6)
    title('Mixed network','FontWeight','normal','FontSize',7);
ax = axes('Position',[0.08,0.14,0.19,0.3]);
labelpanel(ax.Position(1)-0.05,ax.Position(2)+ax.Position(4)-0.01,'e',true);
    plot(f,mean((P_baseline+P_baseline).*(1+0.0*randn(size(P_baseline))),2),'LineWidth',1,'color',grey)
    hold on;
    plot(f,mean((P_baseline+P_baseline_osc).*(1+0.0*randn(size(P_baseline))),2),'LineWidth',1,'color',black)
    plot(f,mean((P_baseline*1.25+0.75*P_crit).*(1+0.0*randn(size(P_baseline))),2),'LineWidth',1,'color',red)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.5,50]);
    xticks([0.5,5,50]);
    ylim(10.^[-16,-14.5]);
    xticklabels([0.5,5,50]);
    xlabel('Frequency (Hz)');
    ylabel('log power')
    % ylabel(['PSD (' char(956) 'V^2/Hz)']);
    yticks([]);
    text(0.6,0.85*2.7*10^-16,'m=0 (ap.)','color',grey,'FontSize',6)
    text(0.6,0.85*2*10^-16,'m=0 (+3 Hz)','color',black,'FontSize',6)
    text(0.6,0.85*1.5*10^-16,'m=0.98 (ap.)','color',red,'FontSize',6)
    title('Freq. res. = 0.1 Hz','FontWeight','normal','FontSize',7);

load('C:\Users\brake\Documents\temp\analzye_simulations_2.mat');
P1(:,:,7:10) = nan;
P_baseline = 1.2*P1(:,:,5);
P_tau = P1(:,:,2);
P_crit = P1(:,:,3);
P_tau_crit = P1(:,:,10);
P_baseline_osc = P1(:,:,6);
P_tau_osc = P1(:,:,8);
P_crit_osc = P1(:,:,9);

asynch = mean((P_baseline+P_baseline).*(1+0.3*randn(size(P_baseline))),2);
osc = mean((P_baseline+P_baseline_osc).*(1+0.3*randn(size(P_baseline))),2);
crit = mean((P_baseline*1.25+0.75*P_crit).*(1+0.3*randn(size(P_baseline))),2);
tau = mean((P_baseline*1.25+0.75*P_tau).*(1+0.3*randn(size(P_baseline))),2);
ax = axes('Position',[0.42,0.14,0.19,0.3]);
labelpanel(ax.Position(1)-0.05,ax.Position(2)+ax.Position(4)-0.01,'f',true);
    plot(f,asynch,'LineWidth',1,'color',grey)
    hold on;
    plot(f,osc,'LineWidth',1,'color',black)
    plot(f,crit,'LineWidth',1,'color',red)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.5,50]);
    xticks([0.5,5,50]);
    ylim(10.^[-16,-14.5]);
    xticklabels([0.5,5,50]);
    xlabel('Frequency (Hz)');
    ylabel('log power')
    % ylabel(['PSD (' char(956) 'V^2/Hz)']);
    yticks([]);
    text(0.6,0.85*2.7*10^-16,'m=0 (ap.)','color',grey,'FontSize',6)
    text(0.6,0.85*2*10^-16,'m=0 (+3 Hz)','color',black,'FontSize',6)
    text(0.6,0.85*1.5*10^-16,'m=0.98 (ap.)','color',red,'FontSize',6)
    title('Freq. res. = 0.5 Hz','FontWeight','normal','FontSize',7);
ax = axes('Position',[0.76,0.14,0.19,0.3]);
labelpanel(ax.Position(1)-0.05,ax.Position(2)+ax.Position(4)-0.01,'g',true);
    plot(f,asynch,'LineWidth',1,'color',grey)
    hold on;
    plot(f,osc,'LineWidth',1,'color',black)
    plot(f,tau,'LineWidth',1,'color',blue)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.5,50]);
    xticks([0.5,5,50]);
    ylim(10.^[-16,-14.5]);
    xticklabels([0.5,5,50]);
    xlabel('Frequency (Hz)');
    ylabel('log power')
    % ylabel(['PSD (' char(956) 'V^2/Hz)']);
    yticks([]);
    text(0.6,0.85*2.7*10^-16,'m=0 (ap.)','color',grey,'FontSize',6)
    text(0.6,0.85*2*10^-16,'m=0 (+3 Hz)','color',black,'FontSize',6)
    text(0.6,0.85*1.5*10^-16,'m=0 (\tau =30)','color',blue,'FontSize',6)
    title('Freq. res. = 0.5 Hz','FontWeight','normal','FontSize',7);

gcaformat(gcf)