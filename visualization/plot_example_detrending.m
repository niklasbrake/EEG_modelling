data = load('E:\Research_Projects\004_Propofol\Experiments\scalp_EEG\analyzed_data\Cz_multitaper_mean.mat');

ptIdx = 1;
temp = data.psd(1,:,ptIdx);
tIdcs = find(and(~isnan(temp),temp~=0));
fIdcs = find(data.freq<300);
ff = data.freq(fIdcs);

pp = data.psd(fIdcs,tIdcs(1),ptIdx);
pp(and(ff>55,ff<65)) = nan;
pp = fillgaps(pp,5);
plotDetrend(ff,pp);

pp = data.psd(fIdcs,tIdcs(end),ptIdx);
pp(and(ff>55,ff<65)) = nan;
pp = fillgaps(pp,5);
plotDetrend(ff,pp);



load('E:\Research_Projects\004_Propofol\Experiments\scalp_EEG\raw_data\timeInformation.mat','timeInfo');
infusionTime = timeInfo.infusion_onset-timeInfo.object_drop;

ptIdx = 12;
tIdcs = find(data.time<infusionTime(ptIdx));
pre = nanmedian(data.psd(fIdcs,tIdcs,ptIdx),2);
temp = data.psd(1,:,ptIdx);
tIdcs = find(and(~isnan(temp),temp~=0));
pp = data.psd(fIdcs,tIdcs(end),ptIdx);
pp(and(ff>55,ff<65)) = nan;
pp = fillgaps(pp,5);

fig = figureNB(5.5,5)
subplot(2,2,1);
    plot(ff,pp(:),'k','LineWidth',1); hold on;
    plot(ff,pre(:),'color','r','LineWidth',1);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    ylim([1e-4,1e4])
    yticks([1e-4,1,1e4]);
 subplot(2,2,2);
    plot(ff,10*log10(pp(:)./pre(:)),'k','LineWidth',1); hold on;
    line([1,150],[0,0],'color','r','linestyle','-','LineWidth',1)
    set(gca,'xscale','log');
    % set(gca,'yscale','log');
    xlim([1,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel(['Power (dB)'])
    ylim([-10,20])
gcaformat(fig)

function plotDetrend(ff,pp)
    [synPre,synFun] = synDetrend(ff,pp,1,'exp2');
    [oofPre,oofFun] = getFOOOF(ff(ff<50),pp(ff<50),false);
    fig = figureNB(5.5,5)
    subplot(2,2,1);
        plot(ff,pp(:),'k','LineWidth',1); hold on;
        plot(ff,pre(:),'k','LineWidth',1); hold on;
        % plot(ff,10.^(synFun(ff,synPre)),'color','r','LineWidth',1);
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        xlim([1,100]);
        xticks([1,10,100]);
        xticklabels([1,10,100]);
        xlabel('Frequency (Hz)');
        ylabel(['PSD (' char(956) 'V^2/Hz)'])
        ylim([1e-4,1e4])
        yticks([1e-4,1,1e4]);
    subplot(2,2,2);
        plot(ff,10*log10(pp(:))-10*synFun(ff,synPre),'k','LineWidth',1); hold on;
        line([1,150],[0,0],'color','r','linestyle','-','LineWidth',1)
        set(gca,'xscale','log');
        % set(gca,'yscale','log');
        xlim([1,100]);
        xticks([1,10,100]);
        xticklabels([1,10,100]);
        xlabel('Frequency (Hz)');
        ylabel(['Power (dB)'])
        ylim([-10,20])
    subplot(2,2,3);
        plot(ff,pp(:),'k','LineWidth',1); hold on;
        plot(ff,oofFun(ff,oofPre),'color','r','LineWidth',1);
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        xlim([1,100]);
        xticks([1,10,100]);
        xticklabels([1,10,100]);
        xlabel('Frequency (Hz)');
        ylabel(['PSD (' char(956) 'V^2/Hz)'])
        ylim([1e-4,1e4])
        yticks([1e-4,1,1e4]);
    subplot(2,2,4);
        plot(ff,10*log10(pp(:))-10*log10(oofFun(ff,oofPre)),'k','LineWidth',1); hold on;
        line([1,150],[0,0],'color','r','linestyle','-','LineWidth',1)
        set(gca,'xscale','log');
        % set(gca,'yscale','log');
        xlim([1,100]);
        xticks([1,10,100]);
        xticklabels([1,10,100]);
        xlabel('Frequency (Hz)');
        ylabel(['Power (dB)'])
        ylim([-10,20])
    gcaformat(fig);
end