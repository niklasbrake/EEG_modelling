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

fig = figureNB(18.3,8);
axes('Position',[0.037,0.55,0.47,0.3]);
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

axes('Position',[0.05,0.11,0.19,0.34]);
    imagesc(time,freq,nanmedian(dB,3));
    axis xy
    ylim([0.5,40])
    CB = colorbar;
    CB.Location = 'eastoutside';
    CB.Position(1) = 0.247;
    CB.Position(3) = 0.015;
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

axes('Position',[0.35,0.11,0.12,0.34]);
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

sim = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\unitary_spectra_propofol_13tau.mat');
axes('Position',[0.865,0.6,0.12,0.32]);
    plot(sim.f,mean(sim.P_baseline_13,2),'LineWidth',1,'color','k')
    hold on;
    plot(sim.f,mean(sim.P_propofol_13,2),'LineWidth',1,'color',red)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    gcaformat;
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    text(0.7,10.^0.11*1e-16,'Baseline','fontsize',7);
    text(7.7,40*1e-16,'Propofol','fontsize',7,'color',red);
    txt = text(0.7,8e-14,'Simulated spectra','fontsize',7);
    xlim([0.5,100])
    ylim([1e-18,1e-13])

axes('Position',[0.55,0.6,0.1,0.22]);
    t = linspace(-20,200,1e3);
    tau1 = 23;
    tau2 = 2;
    fun2 = @(t,tau1) -(exp(-t/tau1) - exp(-t/tau2)).*(t>=0);
    fun = @(t,tau1) fun2(t,tau1)./max(abs(fun2(t,tau1)));
    A = linspace(1,3,5)/3;
    plot(t,A(1)*fun(t,A(1)*tau1),'LineWidth',1.5,'color','k'); hold on;
    plot(t,A(end)*fun(t,A(end)*tau1),'LineWidth',1.5,'color',red); hold on;
    xlim([t(1),t(end)]);
    ylim([-1,0]);
    axis off;
    gam = 0.15;
    text(70,-0.3,'\tau_{decay}','fontsize',7);
    text(90,-0.47,'+200%','fontsize',7,'color',red);
    text(70,-0.73,'g_{max}','fontsize',7);
    text(90,-0.9,'+200%','fontsize',7,'color',red);
    text(90,0.2,'GABA current','FontSize',7,'HorizontalAlignment','center');
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
        n = poissrnd(8);
        ids = [ids;i*ones(n,1)];
        ts = [ts;rand(n,1)];
    end
    hold on;
    R = raster(ids-10,ts,fig);
    R.LineWidth = 1;
    R.Color = red;
    ylim([-9,10]);
    xlim([0,1]);
    axis off;
    text(0,10+0.2*19,'Firing rate','FontSize',7,'HorizontalAlignment','left');
    text(0.63,10+0.2*19,['(' char(8722) '10%)'],'FontSize',7,'HorizontalAlignment','left','color',red);


load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\baseline_preLOC_fits.mat')

[full_model,synFun] = fittingmodel;
for i = 1:14
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
    text(0.6,3e-2,'Baseline (pt.1)','FontSize',7)
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
    text(0.6,3e-2,'Pre-LOC (pt.1)','FontSize',7,'color','k')

axes('Position',[0.865,0.11,0.12,0.34]);
%{
    plotwitherror(freq,pSynPre,'CI','color','k','LineWidth',1);
    plotwitherror(freq,pSynPost,'CI','color',red,'LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xticklabels([1,10,100]);
    xlim([0.5,100])
    xticks([1,10,100])
    ylim([1e-2,1e3])
    text(0.7,10.^0.11,'Baseline','fontsize',7);
    text(7.7,40,'Pre-LOC','fontsize',7,'color',red);
    gcaformat
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
%}



folder = 'E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\rescaled_manual\fitted';
load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\rescaled_Sep2022\rescaled_data.mat')
for i = 1:14
load(fullfile(folder,['pt' int2str(i) '_rescaled_28-Sep-2022.mat']));
params(:,:,i) = pars;
end
tau = squeeze(params(1,:,:))*1e3;

    plotwitherror(linspace(-1.5,0.5,200),squeeze(params(1,:,:))*1e3,'SE')
    gcaformat;



labelpanel(0.0075,0.885,'a',true);
labelpanel(0.0075,0.415,'b',true);
labelpanel(0.3,0.415,'c',true);
labelpanel(0.525,0.885,'d',true);
labelpanel(0.8,0.885,'e',true);
labelpanel(0.5,0.415,'f',true);
labelpanel(0.8,0.415,'g',true);


A = labelpanel(0.545,0.88,'');
A.Position(3) = 0.2;
A.FontSize = 7;
A.FontWeight = 'normal';
A.String = 'Simulated effects of propofol';
A.Color = red;





return;
[full_model,synFun] = fittingmodel;
tRescaled = linspace(-1.5,0.5,200);
folder = 'E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\rescaled_manual\fitted';
params = [];
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
figureNB;
    clrs = clrsPT.sequential(N+m);
    clrs = clrs(1+m:end,:);
    idcs = find(and(tRescaled>-1,tRescaled<0));
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
gcaformat(fig)
axes('Position',[0.6,0.8,0.14,0.03]);
    imagesc(1:N)
    colormap(gca,clrs)
    xticks([1,N])
    yticks([]);
    xticklabels({});
    axis xy;
    xticks([]);
    text(1,0.35,'Infusion','FontSize',6,'HorizontalAlignment','center','VerticalAlignment','top')
    text(N,0.35,'LOC','FontSize',6,'HorizontalAlignment','center','VerticalAlignment','top')
    text((1+N)/2,1.5,'Rescaled time','FontSize',6,'HorizontalAlignment','center','VerticalAlignment','bottom')

    % line(get(gca,'xlim'),1.5*[1,1],'color','k');
    % line((N+0.5)*[1,1],get(gca,'ylim'),,'color','k');



