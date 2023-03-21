eeg_example = load('C:\Users\brake\Documents\GitHub\Propofol2021-private\data\sampleTimeSeries1.mat');

load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\raw\timeInformation.mat')
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\analyzed\Cz_multitaper_mean.mat')

% Compute Baseline spectrum
for i = 1:14
    pre(:,i) = nanmedian(psd(:,and(time>=t0(i)-10,time<t0(i)),i),2);
    post(:,i) = nanmedian(psd(:,and(time>=-10,time<0),i),2);
end

iNoise = find(and(freq>55,freq<65));
pre(iNoise,:) = nan;
post(iNoise,:) = nan;
% Convert to Baseline-normalized decibels
dB = log10(psd);


red = clrsPT.qualitative_CM.red;
clrs = clrsPT.sequential(7);
blue = clrsPT.qualitative_CM.blue;

fig = figureNB(8.9,10);
axes('Position',[0.067,0.76,0.91,0.2]);
    plot(downsample(eeg_example.time,1),downsample(eeg_example.timedomain,1),'LineWidth',0.2,'color',blue); hold on;
    xlim([eeg_example.time(1)-10,eeg_example.time(end)]);
    ylim([-100,100]);
    line([eeg_example.time(1)-10,eeg_example.time(1)-10],[-25,25],'color','k','linewidth',1);
    line([eeg_example.time(1),eeg_example.time(1)+60],[-60,-60],'color','k','linewidth',1);
    text(eeg_example.time(1)-15,0,['50 ' char(956) 'V'],'VerticalAlignment','bottom','HorizontalAlignment','center','fontsize',7,'color','k','Rotation',90);
    text(eeg_example.time(1)+30,-65,'60 s','VerticalAlignment','top','HorizontalAlignment','center','fontsize',7,'color','k');
    scatter(0,90,10,'vk','filled');
    text(0,105,sprintf('LOC'),'FontSize',7,'VerticalAlignment','bottom','HorizontalAlignment','center','color','k');
    line([-200,-1],[75,75],'color',red,'linewidth',2);
    text(-100,80,'Propofol infusion','color',red,'FontSize',7,'VerticalAlignment','bottom','HorizontalAlignment','center');
    axis off;

axes('Position',[0.11,0.465,0.375,0.26]);
    imagesc(time,freq,nanmedian(dB,3));
    axis xy
    ylim([0.5,40])
    CB = colorbar;
    CB.Location = 'eastoutside';
    CB.Position(1) = 0.5;
    CB.Position(3) = 0.03;
    set(gca,'CLim',[-2,3])
    CM = clrsPT.iridescent(50);
    idcs = interp1(linspace(0,1,50),1:50,linspace(-0.3,1.5,50),'nearest','extrap');
    colormap(flip(CM(idcs,:)));
    set(gca,'fontsize',7);
    xlabel('LOC-aligned time (s)');
    ylabel('Frequency (Hz)');
    line([0,0],get(gca,'ylim'),'color','k','linestyle','-','linewidth',1);
    line([0,0],get(gca,'ylim'),'color','w','linestyle','--','linewidth',1);
    xlim([-300,60]);
    CB.Ticks = [-2,0,2,4];
    CB.TickLabels = {'10^{-2}','10^0','10^2','10^{4}'};
    CB.FontSize = 7;

axes('Position',[0.71,0.465,0.25,0.26]);
    plotwitherror(freq,log10(pre),'CI','linewidth',1,'color','k');
    plotwitherror(freq,log10(post),'CI','linewidth',1,'color',red);
    xlim([0.5,300])
    set(gca,'xscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlim([0.5,100]);
    gcaformat;
    xlabel('Frequency (Hz)');
    yl = ylabel(['PSD (' char(956) 'V^2/Hz)']);
    yl.Position(1:2) = [0.1,0.36];
    ylim([-2,3])
    yticks([-2,0,2,4]);
    yticklabels({'10^{-2}','10^0','10^2','10^{4}'});
    text(0.7,0.11,'Baseline','fontsize',7);
    text(2.3,2.6,'Pre-LOC (0-10 s)','fontsize',7,'color',red);


load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\baseline_preLOC_fits.mat')

[full_model,synFun] = fittingmodel;
for i = 1:14
    pSynPre(:,i) = 10.^synFun(freq,synPre(:,i));
    pSynPost(:,i) = 10.^synFun(freq,synPost(:,i));
end

ptIdx = 1;
axes('Position',[0.11,0.09,0.22,0.26]);
    plot(freq,pre(:,ptIdx),'color','k','LineWidth',1);
    hold on;
    plot(freq,10.^synFun(freq,synPre(:,ptIdx)),'color',blue,'linestyle',':','LineWidth',1);
    plot(freq,10.^full_model(freq,synPre(:,ptIdx)),'linestyle','-','linewidth',1,'color',blue);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlim([0.5,100]);
    gcaformat;
    xl = xlabel('Frequency (Hz)');
    xl.Position(1:2) = [150,1.1e-3];
    yl = ylabel(['PSD (' char(956) 'V^2/Hz)'])
    yl.Position(1) = 0.1;
    ylim([1e-2,1e3])
    set(get(gca,'yaxis'),'MinorTick','off');
    text(0.6,3e-2,'Baseline (pt.1)','FontSize',7)
axes('Position',[0.36,0.09,0.22,0.26]);
    plot(freq,post(:,ptIdx),'color','k','LineWidth',1);
    hold on;
    plot(freq,10.^synFun(freq,synPost(:,ptIdx)),'color',blue,'linestyle',':','LineWidth',1);
    plot(freq,10.^full_model(freq,synPost(:,ptIdx)),'linestyle','-','linewidth',1,'color',blue);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlim([0.5,100]);
    gcaformat;
    ylim([1e-2,1e3])
    set(get(gca,'yaxis'),'MinorTick','off');
    yticklabels({});
    text(0.6,3e-2,'Pre-LOC (pt.1)','FontSize',7,'color','k')





folder = 'E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\rescaled_manual\fitted';
load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\rescaled_Sep2022\rescaled_data.mat')
for i = 1:14
load(fullfile(folder,['pt' int2str(i) '_rescaled_28-Sep-2022.mat']));
params(:,:,i) = pars;
end
tau = squeeze(params(1,:,:))*1e3;

axes('Position',[0.71,0.09,0.25,0.26]);
    plotwitherror(linspace(-1.5,0.5,200),squeeze(params(1,:,:))*1e3,'SE','color','k')
    xticks([-1,0]);
    xticklabels({'Infusion','LOC'})
    xlabel('Rescaled time')
    ylabel('\tau_2 (ms)')
    xlim([-1.4,0.4])
    gcaformat


gcaformat(gcf)

labelpanel(0.01,0.95,'a',true);
labelpanel(0.01,0.7,'b',true);
labelpanel(0.61,0.7,'c',true);
labelpanel(0.01,0.325,'d',true);
labelpanel(0.61,0.325,'e',true);;



