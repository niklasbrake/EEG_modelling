dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';
load(fullfile(dataFolder,'data_time_information.mat'));
t0 = timeInfo.infusion_onset-timeInfo.object_drop;
load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\raw\time_series_all_channels.mat')

for i = 1:14
    X = TimeDomainAligned(:,2,i);
    % X = fillgaps(X,1e4); X = X(2e5:6e5);
    X = X(2e5:6e5);
    t = Time(2e5:6e5);
    [freq,time,psd(:,:,i)] = eegfft(t,X,10,9.5);
end



for i = 1:14
    idcs0 = find(and(time<t0(i),time>=-Inf));
    P0(:,i) = nanmedian(psd(:,idcs0,i),2);
    idcs = find(and(time>=-15,time<-5));
    P1(:,i) = nanmedian(psd(:,idcs,i),2);
    idcs2 = find(and(time>=5,time<15));
    P2(:,i) = nanmedian(psd(:,idcs2,i),2);
end
% P2(:,11) = nan;

idcs = find(and(freq>55,freq<65));
P0(idcs,:) = nan; P0 = log10(P0);
P1(idcs,:) = nan; P1 = log10(P1);
P2(idcs,:) = nan; P2 = log10(P2);

figureNB(14,4.5);
subplot(1,3,1);
    plotwitherror(freq,P0,'CI','color','k','LineWidth',1);
    set(gca,'xscale','log')
    % set(gca,'yscale','log')
    xlim([0.1,100])
    ylim([-3,4])
    text(0.15,-1.3,'Baseline','FontSize',7,'color','k');
    % text(0.15,-1.3/4,'pre-LOC (-15 to -5 s)','FontSize',7,'color','b');
    % text(0.15,-1.3/16,'post-LOC (5 to 15 s)','FontSize',7,'color','r');
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)']);;
    xticks([0.1,1,10,100])
    xticklabels([0.1,1,10,100]) 
    yticks([-2:2:4]);
    yticklabels({'10^{-2}','10^0','10^2','10^4'});
    gcaformat
subplot(1,3,2);
    plot(freq,nanmean(P0,2),'color','k','LineWidth',1);
    hold on;
    plotwitherror(freq,P1,'CI','color','b','LineWidth',1);
    set(gca,'xscale','log')
    % set(gca,'yscale','log')
    xlim([0.1,100])
    ylim([-3,4])
    text(0.15,-1.3,'Baseline','FontSize',7,'color','k');
    text(0.15,-1.9,'pre-LOC (-15 to -5 s)','FontSize',7,'color','b');
    % text(0.15,-1.3/16,'post-LOC (5 to 15 s)','FontSize',7,'color','r');
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)']);;
    xticks([0.1,1,10,100])
    xticklabels([0.1,1,10,100]) 
    yticks([-2:2:4]);
    yticklabels({'10^{-2}','10^0','10^2','10^4'});
    gcaformat
subplot(1,3,3);
    plot(freq,nanmean(P0,2),'color','k','LineWidth',1);
    hold on;
    plot(freq,nanmean(P1,2),'color','b','LineWidth',1);
    plotwitherror(freq,P2,'CI','color','r','LineWidth',1);
    set(gca,'xscale','log')
    % set(gca,'yscale','log')
    xlim([0.1,100])
    ylim([-3,4])
    text(0.15,-1.3,'Baseline','FontSize',7,'color','k');
    text(0.15,-1.9,'pre-LOC (-15 to -5 s)','FontSize',7,'color','b');
    text(0.15,-2.5,'post-LOC (5 to 15 s)','FontSize',7,'color','r');
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)']);;
    xticks([0.1,1,10,100])
    xticklabels([0.1,1,10,100]) 
    yticks([-2:2:4]);
    yticklabels({'10^{-2}','10^0','10^2','10^4'});
    gcaformat


figureNB;
for i = 1:14
    subplot(2,7,i);
    plot(freq,P0(:,i))
    hold on;
    plot(freq,P1(:,i))
    plot(freq,P2(:,i))
    set(gca,'xscale','log');
    xlim([0.1,100]);
    ylim([-3,4]);
end