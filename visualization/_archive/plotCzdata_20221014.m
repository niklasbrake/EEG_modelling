eeg_example = load('C:\Users\brake\Documents\GitHub\Propofol2021-private\data\sampleTimeSeries1.mat');

load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\raw\timeInformation.mat')
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\analyzed\Cz_multitaper_mean.mat')

% Compute baseline spectrum
for i = 1:14
    pre(:,i) = nanmedian(psd(:,time<t0(i),i),2);
end

for i = 1:14
    % post(:,i) = nanmedian(psd(:,and(time>=0,time<60),i),2);
    post(:,i) = nanmedian(psd(:,and(time>=-30,time<0),i),2);
end

iNoise = find(and(freq>55,freq<65));
pre(iNoise,:) = nan;
post(iNoise,:) = nan;
% Convert to baseline-normalized decibels
dB = log10(psd);


red = clrsPT.qualitative_CM.red;
clrs = clrsPT.sequential(7);
% blue = [0.0667,0.2,0.8];
blue = clrsPT.qualitative_CM.blue;
fig = figureNB(18.3,8);
axes('Position',[0.037,0.55,0.47,0.3]);
% subplot(2,2,1);
    plot(downsample(eeg_example.time,1),downsample(eeg_example.timedomain,1),'LineWidth',0.2,'color',blue); hold on;
    xlim([eeg_example.time(1)-10,eeg_example.time(end)]);
    ylim([-100,100]);
    line([eeg_example.time(1)-10,eeg_example.time(1)-10],[-25,25],'color','k','linewidth',1);
    line([eeg_example.time(1),eeg_example.time(1)+60],[-60,-60],'color','k','linewidth',1);
    text(eeg_example.time(1)-15,0,['50 ' char(956) 'V'],'VerticalAlignment','bottom','HorizontalAlignment','center','fontsize',7,'color','k','Rotation',90);
    text(eeg_example.time(1)+30,-65,'60 s','VerticalAlignment','top','HorizontalAlignment','center','fontsize',7,'color','k');
    scatter(0,90,10,'vk','filled'); text(7.5,105,sprintf('LOC'),'FontSize',7,'VerticalAlignment','bottom','HorizontalAlignment','center','color','k');
    line([-200,-1],[75,75],'color',red,'linewidth',2);
    text(-100,80,'Propofol Infusion','color',red,'FontSize',7,'VerticalAlignment','bottom','HorizontalAlignment','center');
    axis off;
axes('Position',[0.05,0.11,0.19,0.34]);
% subplot(2,4,5);
    imagesc(time,freq,nanmedian(dB,3));
    axis xy
    ylim([0.5,40])
    CB = colorbar;
    CB.Location = 'eastoutside';
    CB.Position(1) = 0.247;
    CB.Position(3) = 0.015;
    % CB.Position(3) = 0.1;
    % CB.Title.String = '(dB)';
    % set(gca,'CLim',[-10,30])
    % set(gca,'CLim',[0,1])
    CM = jet(1e3);
    CM = CM(1:floor(1e3/6*5),:);
    set(gca,'CLim',[-2,3])
    colormap(CM)
    set(gca,'fontsize',7);
    xlabel('LOC-aligned time (s)');
    ylabel('Frequency (Hz)');
    line([0,0],get(gca,'ylim'),'color','k','linestyle','-','linewidth',1);
    line([0,0],get(gca,'ylim'),'color','w','linestyle','--','linewidth',1);
    xlim([-300,60]);
    CB.Ticks = [-2,0,2,4];
    CB.TickLabels = {'10^{-2}','10^0','10^2','10^{4}'};
    CB.FontSize = 7;
axes('Position',[0.35,0.11,0.12,0.34]);
% subplot(2,4,6);
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
    text(0.7,0.11,'baseline','fontsize',7);
    text(2.3,2.6,'pre-LOC (0-30 s)','fontsize',7,'color',red);



propo = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\expected_PSD_propofol.mat');
baseline = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\expected_PSD_baseline.mat');

axes('Position',[0.865,0.6,0.12,0.32]);
% subplot(2,4,4);
    plotwitherror(baseline.f,baseline.p,'CI','color','k','LineWidth',1);
    plotwitherror(propo.f,1.75*propo.p,'CI','color',red,'LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlim([0.5,100]);
    gcaformat;
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    text(0.7,9e-19,'baseline','fontsize',7);
    text(0.8,1.9e-17,'Propofol simulation','fontsize',7,'color',red);

axes('Position',[0.55,0.6,0.1,0.22]);
    t = linspace(-20,200,1e3);
    tau1 = 23;
    tau2 = 2;
    fun2 = @(t,tau1) -(exp(-t/tau1) - exp(-t/tau2)).*(t>=0);
    fun = @(t,tau1) fun2(t,tau1)./max(abs(fun2(t,tau1)));
    A = linspace(1,3,5);
    for i = 1:5
        plot(t,A(i)*fun(t,A(i)*tau1),'LineWidth',1.5,'color',clrs(i,:)); hold on;
    end
    xlim([t(1),t(end)]);
    ylim([-3,0]);
    axis off;
    text(-10,0.2,'\tau \rightarrow','fontsize',7)
    text(60,0.2,'2.5\tau','fontsize',7,'color',red)

axes('Position',[0.67,0.6,0.1,0.22]);
    ids = []; ts = [];
    for i = 1:10
        n = poissrnd(15);
        ids = [ids;i*ones(n,1)];
        ts = [ts;rand(n,1)];
    end
    R = raster(ids,ts,fig);
    R.LineWidth = 1;
    ids = []; ts = [];
    for i = 1:10
        n = poissrnd(6);
        ids = [ids;i*ones(n,1)];
        ts = [ts;rand(n,1)];
    end
    hold on;
    R = raster(ids-10,ts,fig);
    R.LineWidth = 1;
    R.Color = red;
    ylim([-9,10]);
    axis off;
    text(0,14,'\lambda \rightarrow','fontsize',7)
    text(0.3,14,'0.5\lambda','fontsize',7,'color',red)

load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\pre_post_fits.mat')
pSynPre=[];
pSynPost=[];
[full_model,synFun] = fittingmodel;
for i =1:14
    pSynPre(:,i) = 10.^synFun(freq,synPre(:,i));
    pSynPost(:,i) = 10.^synFun(freq,synPost(:,i));
end

ptIdx = 1;
axes('Position',[0.55,0.11,0.12,0.34]);
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
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    ylim([1e-2,1e3])
    set(get(gca,'yaxis'),'MinorTick','off');
    text(0.6,3e-2,'baseline (pt.1)','FontSize',7)
axes('Position',[0.68,0.11,0.12,0.34]);
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
    text(0.6,3e-2,'pre-LOC (pt.1)','FontSize',7,'color','k')
axes('Position',[0.865,0.11,0.12,0.34]);
    plotwitherror(freq,pSynPre,'SE','color','k','LineWidth',1); hold on;
    plotwitherror(freq,pSynPost,'SE','color',red,'LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlim([0.5,100]);
    gcaformat;
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    ylim([1e-2,1e3])
    set(get(gca,'yaxis'),'MinorTick','off');
    text(1,1,'baseline','FontSize',7);
    text(2,82,'pre-LOC','FontSize',7,'color',red);
gcaformat(fig)

CM = clrsPT.iridescent(50);
idcs = interp1(linspace(0,1,50),1:50,linspace(-0.3,1.5,50),'nearest','extrap');
colormap(flip(CM(idcs,:)));

labelpanel(0,0.885,'a',true);
labelpanel(0,0.415,'b',true);
labelpanel(0.3,0.415,'c',true);
labelpanel(0.53,0.885,'d',true);
labelpanel(0.8,0.885,'e',true);
labelpanel(0.5,0.415,'f',true);
labelpanel(0.8,0.415,'g',true);


A = labelpanel(0.55,0.88,'');
A.Position(3) = 0.2;
A.FontSize = 7;
A.FontWeight = 'normal';
A.String = 'Simulated effects of propofol';
A.Color = red;






















load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\unitary_spectra_propofol.mat');


[full_model,synFun] = fittingmodel;
tRescaled = linspace(-1.5,0.5,200);
folder = 'E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\rescaled_manual\fitted';
for i = 1:14
    load(fullfile(folder,['pt' int2str(i) '_rescaled_28-Sep-2022.mat']));
    params(:,:,i) = pars;
end
ff = 0.5:0.5:150;
for i = 1:14
    for j = 1:size(params,2)
        ptBL(:,j,i) = synFun(ff,params(:,j,i));
    end
end
N = 7;
m = 5;

clrs = clrsPT.sequential(length(P)+m);
clrs = clrs(m:end,:);
figureNB(5,11);
subplot(2,1,1);
    for i = 1:length(P)
        plot(f,mean(P{i},2),'LineWidth',1,'color',clrs(i,:))
        hold on;
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    gcaformat;
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlim([0.5,100])
    ylim([1e-18,1e-13])

subplot(2,1,2);
    clrs = clrsPT.sequential(N+m);
    clrs = clrs(1+m:end,:);

    idcs = find(and(tRescaled>-1,tRescaled<.5));
    k = floor(length(idcs)/N);
    for i = 1:N
        BL = mean(mean(ptBL(:,idcs(k*(i-1)+1:k*i),:),2),3);
        plot(ff,10.^BL,'color',clrs(i,:),'LineWidth',1); hold on;
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.5,100])
    xticks([1,10,100])
    ylim([1e-2,1e3])
    gcaformat
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    colormap(clrs)
    CB = colorbar('location','south');
    CB.TickLabels = {'Infusion','LOC'};
    CB.Label.String = 'Rescaled time';
    CB.TickDirection = 'out';
    CB.TickLength = 0.05;
    CB.Ticks = [0.5/N,1-0.5/N];