load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\pre.mat');
load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\fits.mat');

ptIdx = 1;
preExample = 10.^nanmedian(rescaled.psd(:,rescaled.t<-1,ptIdx),2);
postExample = 10.^nanmedian(rescaled.psd(:,rescaled.t>0,ptIdx),2);
[oofPre,oofFun] = getFOOOF(freq(freq<50),preExample(freq<50),false);
[oofPost,oofFun] = getFOOOF(freq(freq<50),postExample(freq<50),false);
[synPre,synFun] = synDetrend(freq(freq<100),preExample(freq<100),3,'exp2');
[synPost,synFun,full_model] = synDetrend(freq(freq<100),postExample(freq<100),3,'exp2');


for i = 1:14
    oofFit = []; synFit = [];
    for j = 1:size(aligned.psd,2)
        oofFit(:,j) = oofFun(freq,aligned.oofPars(:,j,i));
        synFit(:,j) = synFun(freq,aligned.synPars(1:4,j,i));
    end
    aligned.syn(:,:,i) = aligned.psd(:,:,i)-synFit;
    aligned.oof(:,:,i) = aligned.psd(:,:,i)-log10(oofFit);
    aligned.pre(:,:,i) = aligned.psd(:,:,i)-log10(pre(:,i));

    oofFit = []; synFit = [];
    for j = 1:size(rescaled.psd,2)
        oofFit(:,j) = oofFun(freq,rescaled.oofPars(:,j,i));
        synFit(:,j) = synFun(freq,rescaled.synPars(1:4,j,i));
        fitted(:,j) = synFun(freq,rescaled.synPars(1:4,j,i));
    end
    rescaled.syn(:,:,i) = rescaled.psd(:,:,i)-synFit;
    rescaled.fit(:,:,i) = rescaled.psd(:,:,i)-synFit;
    rescaled.oof(:,:,i) = rescaled.psd(:,:,i)-log10(oofFit);
    rescaled.pre(:,:,i) = rescaled.psd(:,:,i)-log10(pre(:,i));

    rescaled.oof(:,:,i) = rescaled.oof(:,:,i)-nanmedian(rescaled.oof(:,rescaled.t<-1,i),2);
    rescaled.syn(:,:,i) = rescaled.syn(:,:,i)-nanmedian(rescaled.syn(:,rescaled.t<-1,i),2);
end


deltaIdx = find(and(freq>=0.5,freq<4));
alphaIdx = find(and(freq>=8,freq<15));
betaIdx = find(and(freq>=15,freq<30));

rescaled.alpha_pre = 10*squeeze(nanmean(rescaled.pre(alphaIdx,:,:)));
rescaled.beta_pre = 10*squeeze(nanmean(rescaled.pre(betaIdx,:,:)));
rescaled.delta_pre = 10*squeeze(nanmean(rescaled.pre(deltaIdx,:,:)));

rescaled.alpha_oof = 10*squeeze(nanmean(rescaled.oof(alphaIdx,:,:)));
rescaled.beta_oof = 10*squeeze(nanmean(rescaled.oof(betaIdx,:,:)));
rescaled.delta_oof = 10*squeeze(nanmean(rescaled.oof(deltaIdx,:,:)));



load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\rescaled_Sep2022\rescaled_data.mat')
folder = 'E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\rescaled_manual\fitted';
load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\rescaled_Sep2022\rescaled_data.mat')
for i = 1:14
    load(fullfile(folder,['pt' int2str(i) '_rescaled_28-Sep-2022.mat']));
    params(:,:,i) = pars;
end
ptBL = nan(size(pRescaled));
fitted = [];
for i = 1:14
    for j = 1:size(params,2)
        ptBL(:,j,i) = synFun(ff,params(:,j,i));
        fitted(:,j,i) = full_model(ff,params(:,j,i));
    end
end
rescaled.fit = fitted;
rescaled.syn = pRescaled-ptBL;
rescaled.syn = rescaled.syn-nanmedian(rescaled.syn(:,tRescaled<-1,:),2);
rescaled.alpha_syn = 10*squeeze(nanmean(rescaled.syn(alphaIdx,:,:)));
rescaled.beta_syn = 10*squeeze(nanmean(rescaled.syn(betaIdx,:,:)));
rescaled.delta_syn = 10*squeeze(nanmean(rescaled.syn(deltaIdx,:,:)));

blue = clrsPT.qualitative_CM.blue;
clrs = clrsPT.lines(3);
clrs(2:3,:) = flip(clrs(2:3,:));
fig = figureNB(12,12);
axes('Position',[0.09,0.72,0.11,0.14]);
    plot(freq,preExample,'LineWidth',1,'color',blue); hold on;
    plot(freq,postExample,'k','LineWidth',1); hold on;
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlabel('Frequency (Hz)');
    ylim([1e-2,1e4]);
    yticks(10.^[-2:2:4]);
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlim([0.5,100]);
    xticks([0.5,5,50]);
    xticklabels([0.5,5,50]);
    set(get(gca,'xaxis'),'MinorTick','off');
    set(get(gca,'yaxis'),'MinorTick','off');
axes('Position',[0.29,0.72,0.11,0.14]);
    plot(freq,10*log10(postExample./preExample),'k','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.5,100]);
    line(get(gca,'xlim'),[0,0],'LineWidth',1,'color',blue,'linestyle','-');
    xticks([0.5,5,50]);
    xticklabels([0.5,5,50]);
    set(get(gca,'xaxis'),'MinorTick','off');
    set(get(gca,'yaxis'),'MinorTick','off');
    xlabel('Frequency (Hz)');
    ylim([-10,20]);
    ylabel('Power (dB)');
axes('Position',[0.51,0.72,0.2,0.18]);
    imagesc(aligned.t,freq,10*nanmedian(aligned.pre,3));
    ylim([0.5,50])
    colormap('jet')
    % colorbar;
    axis xy
    xlim([-180,60]);
    xticks([-180:60:60]);
    xticklabels(-3:1)
    ylabel('Frequency (Hz)')
    xlabel('Time rel. LOC (min)')
    set(gca,'CLim',[-5,10])
axes('Position',[0.79,0.72,0.2,0.18]);
    h(1) = plotwitherror(rescaled.t,smoothdata(rescaled.alpha_pre,'movmedian',10),'SE','LineWidth',1,'color',clrs(1,:));
    h(2) = plotwitherror(rescaled.t,smoothdata(rescaled.beta_pre,'movmedian',10),'SE','LineWidth',1,'color',clrs(2,:));
    h(3) = plotwitherror(rescaled.t,smoothdata(rescaled.delta_pre,'movmedian',10),'SE','LineWidth',1,'color',clrs(3,:));
    xlim([-1.25,0.25]);
    ylim([-2.5,15])
    ylabel('Power (dB)');
    xlabel('Rescaled time')
    xticks([-1,0]);
    % xticklabels([0,0.5,1]);
    xticklabels({'Infusion','LOC'})
    line([-1.5,0.5],[0,0],'color','k');
    text(-1.15,13,'Alpha (8-15 Hz)','FontSize',7,'Color',clrs(1,:));
    text(-1.15,11,'Beta (15-30 Hz)','FontSize',7,'Color',clrs(2,:));
    text(-1.15,15,'Delta (0.5-4 Hz)','FontSize',7,'Color',clrs(3,:));

axes('Position',[0.09,0.4,0.11,0.14]);
    plot(freq,oofFun(freq,oofPost),'LineWidth',1,'color',blue);  hold on
    plot(freq,postExample,'k','LineWidth',1);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xticks([0.5,5,50]);
    xticklabels([0.5,5,50]);
    set(get(gca,'xaxis'),'MinorTick','off');
    set(get(gca,'yaxis'),'MinorTick','off');
    xlabel('Frequency (Hz)');
    ylim([1e-2,1e4]);
    yticks(10.^[-2:2:4]);
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
axes('Position',[0.29,0.4,0.11,0.14]);
    plot(freq,10*log10(postExample./oofFun(freq,oofPost)),'k','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.5,100]);
    line(get(gca,'xlim'),[0,0],'LineWidth',1,'color',blue,'linestyle','-');
    xlim([0.5,100]);
    xticks([0.5,5,50]);
    xticklabels([0.5,5,50]);
    set(get(gca,'xaxis'),'MinorTick','off');
    set(get(gca,'yaxis'),'MinorTick','off');
    xlabel('Frequency (Hz)');
    ylim([-10,20]);
    ylabel('Power (dB)');
axes('Position',[0.51,0.4,0.2,0.18]);
    imagesc(aligned.t,freq,10*nanmedian(aligned.oof,3));
    ylim([0.5,50])
    colormap('jet')
    % colorbar;
    axis xy
    xlim([-180,60]);
    xticks([-180:60:60]);
    xticklabels(-3:1)
    ylabel('Frequency (Hz)')
    xlabel('Time rel. LOC (min)')
    set(gca,'CLim',[-5,10])
axes('Position',[0.79,0.4,0.2,0.18]);
    plotwitherror(rescaled.t,smoothdata(rescaled.alpha_oof,'movmedian',10),'SE','LineWidth',1,'color',clrs(1,:));
    plotwitherror(rescaled.t,smoothdata(rescaled.beta_oof,'movmedian',10),'SE','LineWidth',1,'color',clrs(2,:));
    plotwitherror(rescaled.t,smoothdata(rescaled.delta_oof,'movmedian',10),'SE','LineWidth',1,'color',clrs(3,:));
     xlim([-1.25,0.25]);
     ylim([-2.5,6])
     ylabel('Power (dB)');
     xlabel('Rescaled time')
     xticks([-1,0]);
     % xticklabels([0,0.5,1]);
     xticklabels({'Infusion','LOC'})
     line([-1.5,0.5],[0,0],'color','k');
synPost = [0.047,0.0075,-14.25,4.3];
axes('Position',[0.09,0.08,0.11,0.14]);
    plot(freq,postExample,'k','LineWidth',1); hold on;
    plot(freq,10.^synFun(freq,synPost),'LineWidth',1,'color',blue);  hold on
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xticks([0.5,5,50]);
    xticklabels([0.5,5,50]);
    set(get(gca,'xaxis'),'MinorTick','off');
    set(get(gca,'yaxis'),'MinorTick','off');
    xlabel('Frequency (Hz)');
    ylim([1e-2,1e4]);
    yticks(10.^[-2:2:4]);
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
axes('Position',[0.29,0.08,0.11,0.14]);
    plot(freq,10*(log10(postExample)-synFun(freq,synPost)),'k','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.5,100]);
    line(get(gca,'xlim'),[0,0],'LineWidth',1,'color',blue,'linestyle','-');
    xlim([0.5,100]);
    xticks([0.5,5,50]);
    xticklabels([0.5,5,50]);
    set(get(gca,'xaxis'),'MinorTick','off');
    set(get(gca,'yaxis'),'MinorTick','off');
    xlabel('Frequency (Hz)');
    ylim([-10,20]);
    ylabel('Power (dB)');
axes('Position',[0.51,0.08,0.2,0.18]);
    imagesc(aligned.t,freq,10*nanmedian(aligned.syn,3));
    ylim([0.5,50])
    colormap('jet')
    % colorbar;
    axis xy
    xlim([-180,60]);
    xticks([-180:60:60]);
    xticklabels(-3:1)
    ylabel('Frequency (Hz)')
    xlabel('Time rel. LOC (min)')
    set(gca,'CLim',[-5,10])
axes('Position',[0.79,0.08,0.2,0.18]);
    plotwitherror(tRescaled,smoothdata(rescaled.alpha_syn,'movmedian',3),'SE','LineWidth',1,'color',clrs(1,:));
    plotwitherror(tRescaled,smoothdata(rescaled.beta_syn,'movmedian',3),'SE','LineWidth',1,'color',clrs(2,:));
    plotwitherror(tRescaled,smoothdata(rescaled.delta_syn,'movmedian',3),'SE','LineWidth',1,'color',clrs(3,:));
    xlim([-1.25,0.25]);
    ylabel('Power (dB)');
    xlabel('Rescaled time')
    xticks([-1,0]);
    % xticklabels([0,0.5,1]);
    xticklabels({'Infusion','LOC'})
    line([-1.5,0.5],[0,0],'color','k');
    ylim([-2.5,5])

gcaformat(fig);


idcs = [interp1(linspace(0,1,200),1:200,linspace(0,1,50),'nearest','extrap'),interp1(linspace(0,1,200),201:400,linspace(0,1,100),'nearest','extrap')];
CM = clrsPT.iridescent(400);
colormap(flip(CM(idcs,:)));

labelpanel(0.01,0.86,'a',true);
labelpanel(0.44,0.86,'b',true);
labelpanel(0.725,0.86,'c',true);

labelpanel(0.01,0.54,'d',true);
labelpanel(0.44,0.54,'e',true);
labelpanel(0.725,0.54,'f',true);

labelpanel(0.01,0.22,'g',true);
labelpanel(0.44,0.22,'h',true);
labelpanel(0.725,0.22,'i',true);

A = labelpanel(0.05,0.857,'');
A.String = 'Baseline normalized';
A.FontSize = 7;
A.FontWeight = 'normal';
A.Position(3) = 0.3;

A = labelpanel(0.05,0.537,'');
A.String = '1/f detrended';
A.FontSize = 7;
A.FontWeight = 'normal';
A.Position(3) = 0.3;

A = labelpanel(0.05,0.217,'');
A.String = 'Lorentzian detrended';
A.FontSize = 7;
A.FontWeight = 'normal';
A.Position(3) = 0.3;