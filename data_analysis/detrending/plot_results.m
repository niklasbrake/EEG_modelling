dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';

load(fullfile(dataFolder,'data_time_information.mat'));
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

psd = [];
aligned = load(fullfile(dataFolder,'data_aligned_detrended.mat'));
rescaled = load(fullfile(dataFolder,'data_rescaled_detrended.mat'));
freq = rescaled.freq;
%{

ptIdx = 1;
preExample = 10.^nanmedian(rescaled.psd(:,rescaled.time<-1,ptIdx),2);
postExample = 10.^nanmedian(rescaled.psd(:,rescaled.time>0,ptIdx),2);
[oofPre,oofFun] = getFOOOF(freq(freq<50),preExample(freq<50),false);
[oofPost,oofFun] = getFOOOF(freq(freq<50),postExample(freq<50),false);
[synPre,synFun] = synDetrend(freq(freq<100),preExample(freq<100),3,'exp2');
[synPost,synFun,full_model] = synDetrend(freq(freq<100),postExample(freq<100),3,'exp2');

%}


for i = 1:14
    oofFit = []; synFit = [];
    for j = 1:size(rescaled.psd,2)
        oofFit(:,j) = log10(oofFun(freq,rescaled.oofPars(:,j,i)));
        synFit(:,j) = synFun(freq,rescaled.synPars(1:4,j,i));
    end
    rescaled.syn_detrended(:,:,i) = rescaled.psd(:,:,i)-synFit;
    rescaled.oof_detrended(:,:,i) = rescaled.psd(:,:,i)-oofFit;

    rescaled.pre(:,:,i) = rescaled.psd(:,:,i)-nanmedian(rescaled.psd(:,rescaled.time<-1,i),2);
    rescaled.oof(:,:,i) = rescaled.oof_detrended(:,:,i)-nanmedian(rescaled.oof_detrended(:,rescaled.time<-1,i),2);
    rescaled.syn(:,:,i) = rescaled.syn_detrended(:,:,i)-nanmedian(rescaled.syn_detrended(:,rescaled.time<-1,i),2);

    aligned.syn(:,:,i) = aligned.syn(:,:,i)-nanmedian(aligned.syn(:,aligned.time<t0(i),i),2);
end



deltaIdx = find(and(freq>=1,freq<4));
alphaIdx = find(and(freq>=8,freq<15));
betaIdx = find(and(freq>=15,freq<30));
% betaIdx = find(and(freq>=30,freq<55));

rescaled.alpha_pre = 10*squeeze(nanmean(rescaled.pre(alphaIdx,:,:)));
rescaled.beta_pre = 10*squeeze(nanmean(rescaled.pre(betaIdx,:,:)));
rescaled.delta_pre = 10*squeeze(nanmean(rescaled.pre(deltaIdx,:,:)));

rescaled.alpha_oof = 10*squeeze(nanmean(rescaled.oof(alphaIdx,:,:)));
rescaled.beta_oof = 10*squeeze(nanmean(rescaled.oof(betaIdx,:,:)));
rescaled.delta_oof = 10*squeeze(nanmean(rescaled.oof(deltaIdx,:,:)));

rescaled.alpha_syn = 10*squeeze(nanmean(rescaled.syn(alphaIdx,:,:)));
rescaled.beta_syn = 10*squeeze(nanmean(rescaled.syn(betaIdx,:,:)));
rescaled.delta_syn = 10*squeeze(nanmean(rescaled.syn(deltaIdx,:,:)));

blue = clrsPT.qualitative_CM.blue;
clrs = clrsPT.lines(3);
clrs(2:3,:) = flip(clrs(2:3,:));
fig = figureNB(18,8);
axes('Position',[0.07,0.62,0.08,0.25]);
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
axes('Position',[0.2,0.62,0.08,0.25]);
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
    % ylabel('Power (dB)');
    ylabel('Power (dB)')
axes('Position',[0.34,0.62,0.19,0.3]);
    imagesc(aligned.time,freq,10*nanmean(aligned.pre,3));
    ylim([0.5,50])
    CB = colorbar('location','eastoutside');
    CB.Label.String = 'Power (dB)';
    CB.Position(1) = 0.54;
    % colorbar;
    axis xy
    xlim([-180,60]);
    xticks([-180:60:60]);
    xticklabels(-3:1)
    ylabel('Frequency (Hz)')
    xlabel('Time rel. LOC (min)')
    set(gca,'CLim',[-5,10])
ax1 = axes('Position',[0.82,0.62,0.15,0.3]);
    h(1) = plotwitherror(rescaled.time,smoothdata(rescaled.alpha_pre,'movmedian',3),'SE','LineWidth',1,'color',clrs(1,:));
    h(2) = plotwitherror(rescaled.time,smoothdata(rescaled.beta_pre,'movmedian',3),'SE','LineWidth',1,'color',clrs(2,:));
    h(3) = plotwitherror(rescaled.time,smoothdata(rescaled.delta_pre,'movmedian',3),'SE','LineWidth',1,'color',clrs(3,:));
    xlim([-1.25,0.25]);
    ylim([-2.5,15])
    % ylabel('Power (dB)');
    ylabel('Power (dB)')
    xlabel('Rescaled time')
    xticks([-1,0]);
    % xticklabels([0,0.5,1]);
    xticklabels({'Infusion','LOC'})
    line([-1.5,0.5],[0,0],'color','k');
    text(-1.3,9.5,'\alpha (8-15 Hz)','FontSize',6,'Color',clrs(1,:));
    text(-1.3,7,'\beta (15-30 Hz)','FontSize',6,'Color',clrs(2,:));
    text(-1.3,12,'\delta (1-4 Hz)','FontSize',6,'Color',clrs(3,:));

synPost = [0.047,0.0075,-14.25,4.3];
axes('Position',[0.07,0.15,0.08,0.25]);
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
axes('Position',[0.2,0.15,0.08,0.25]);
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
    % ylabel('Power (dB)');
    ylabel('Power (dB)')
axes('Position',[0.34,0.15,0.19,0.3]);
    imagesc(aligned.time,freq,10*nanmean(aligned.syn,3));
    ylim([0.5,50])
    CB = colorbar('location','eastoutside');
    CB.Label.String = 'Power (dB)';
    CB.Position(1) = 0.54;
    % colorbar;
    axis xy
    xlim([-180,60]);
    xticks([-180:60:60]);
    xticklabels(-3:1)
    ylabel('Frequency (Hz)')
    xlabel('Time rel. LOC (min)')
    set(gca,'CLim',[-3,6])
ax2 = axes('Position',[0.82,0.15,0.15,0.3]);
    plotwitherror(rescaled.time,smoothdata(rescaled.alpha_syn,'movmedian',3),'SE','LineWidth',1,'color',clrs(1,:));
    plotwitherror(rescaled.time,smoothdata(rescaled.beta_syn,'movmedian',3),'SE','LineWidth',1,'color',clrs(2,:));
    plotwitherror(rescaled.time,smoothdata(rescaled.delta_syn,'movmedian',3),'SE','LineWidth',1,'color',clrs(3,:));
    xlim([-1.25,0.25]);
    % ylabel('Power (dB)');
    ylabel('Power (dB)')
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




k = 2;
A = 10*squeeze(nanmean(rescaled.pre(alphaIdx,:,:)));
B = 10*squeeze(nanmean(rescaled.pre(betaIdx,:,:)));
D = 10*squeeze(nanmean(rescaled.pre(deltaIdx,:,:)));
rescaled.p = [];
for i = 1:floor(200/k)
    rescaled.p(i,1) = signtest(nanmedian(A((i-1)*k+1:i*k,:)),0,'tail','right');
    rescaled.p(i,2) = signtest(nanmedian(B((i-1)*k+1:i*k,:)),0,'tail','right');
    rescaled.p(i,3) = signtest(nanmedian(D((i-1)*k+1:i*k,:)),0,'tail','right');
end

A = 10*squeeze(nanmean(rescaled.syn(alphaIdx,:,:)));
B = 10*squeeze(nanmean(rescaled.syn(betaIdx,:,:)));
D = 10*squeeze(nanmean(rescaled.syn(deltaIdx,:,:)));
rescaled.p_syn = [];
rescaled.time2 = [];
for i = 1:floor(200/k)
    rescaled.p_syn(i,1) = signtest(nanmedian(A((i-1)*k+1:i*k,:)),0,'tail','right');
    rescaled.p_syn(i,2) = signtest(nanmedian(B((i-1)*k+1:i*k,:)),0,'tail','right');
    rescaled.p_syn(i,3) = signtest(nanmedian(D((i-1)*k+1:i*k,:)),0,'tail','right');
    rescaled.time2(i) = rescaled.time(k*(i-1)+1);
end

dt = rescaled.time(2)-rescaled.time(1);

idx = [2,1,3];
axes(ax1)
    for j = 1:3
        for i = 1:length(rescaled.time2)
            t1 =rescaled.time2(i);
            fill([t1,t1+k*dt,t1+k*dt,t1],[0,0,2,2]*(rescaled.p(i,j)<0.05)+12+2*idx(j),clrs(j,:),'LineStyle','none');
            hold on
        end
    end
    xlim([-1.4,0.4]);
    % ylim([-2.5,18])
    ylim([-8,20])
    % line([0,0],[0,3],'color','k','LineWidth',1);
    xticks([-1,-0.5,0]);
    xticklabels({'Infusion','','LOC'});
    text(0.43,18.2,'*')
    text(0.43,16.2,'*')
    text(0.43,14.2,'*')
axes(ax2);
    for j = 1:3
        for i = 1:length(rescaled.time2)
            t1 =rescaled.time2(i);
            fill([t1,t1+k*dt,t1+k*dt,t1],[0,0,1,1]*(rescaled.p_syn(i,j)<0.05)+6+idx(j),clrs(j,:),'LineStyle','none');
            hold on
        end
    end
    xlim([-1.4,0.4]);
    ylim([-4,10])
    % line([0,0],[0,3],'color','k','LineWidth',1);
    xticks([-1,-0.5,0]);
    xticklabels({'Infusion','','LOC'});
    text(0.43,9.1,'*')
    text(0.43,8.1,'*')
    text(0.43,7.1,'*')



for i = 1:14
    tt = -t0(1)*rescaled.time;
    idcs = find(and(tt>=-10,tt<=0));
    post(:,i) = 10*nanmedian(rescaled.pre(:,idcs,i),2);
end
A = nanmean(post(alphaIdx,:)); signtest(A,0,'tail','right')
B = nanmean(post(betaIdx,:)); signtest(B,0,'tail','right')
D = nanmean(post(deltaIdx,:)); signtest(D,0,'tail','right')
axes('Position',[0.65,0.62,0.12,0.3]);
    plotwitherror(rescaled.freq,post,'CI','color','k','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.5,100]);
    line(get(gca,'xlim'),[0,0],'LineWidth',1,'color','k','linestyle','--');
    xticks([0.5,5,50]);
    xticklabels([0.5,5,50]);
    set(get(gca,'xaxis'),'MinorTick','off');
    set(get(gca,'yaxis'),'MinorTick','off');
    xlabel('Frequency (Hz)');
    % ylim([-10,20]);
    % ylabel('Power (dB)');
    ylabel('Power (dB)')
    ylim([-8,20]);
    gcaformat;

    fill([8,15,15,8],[18.5,18.5,19,19]-1,clrs(1,:),'LineStyle','none');
    % text(10.^((log10(15)+log10(8))/2),19.5,'*','FontSize',7,'HorizontalAlignment','center');
    fill([15,30,30,15],[18.5,18.5,19,19]-1,clrs(2,:),'LineStyle','none');
    % text(10.^((log10(15)+log10(30))/2),19.5,'*','FontSize',7,'HorizontalAlignment','center');
    fill([1,4,4,1],[18.5,18.5,19,19]-1,clrs(3,:),'LineStyle','none');
    % text(10.^((log10(0.5)+log10(4))/2),19.5,'*','FontSize',7,'HorizontalAlignment','center');

    text(10.^((log10(15)+log10(8))/2),20,'\alpha*','FontSize',7,'HorizontalAlignment','center','color',clrs(1,:));
    text(10.^((log10(15)+log10(30))/2),20,'\beta*','FontSize',7,'HorizontalAlignment','center','color',clrs(2,:));
    text(10.^((log10(1)+log10(4))/2),20,'\delta*','FontSize',7,'HorizontalAlignment','center','color',clrs(3,:));

for i = 1:14
    tt = -t0(1)*rescaled.time;
    idcs = find(and(tt>=-10,tt<=0));
    post_syn(:,i) = 10*nanmedian(rescaled.syn(:,idcs,i),2);
end
A_syn = nanmean(post_syn(alphaIdx,:)); signtest(A_syn,0,'tail','right')
B_syn = nanmean(post_syn(betaIdx,:)); signtest(B_syn,0,'tail','right')
D_syn = nanmean(post_syn(deltaIdx,:)); signtest(D_syn,0,'tail','right')
axes('Position',[0.65,0.15,0.12,0.3]);
    plotwitherror(rescaled.freq,post_syn,'CI','color','k','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.5,100]);
    line(get(gca,'xlim'),[0,0],'LineWidth',1,'color','k','linestyle','--');
    xticks([0.5,5,50]);
    xticklabels([0.5,5,50]);
    set(get(gca,'xaxis'),'MinorTick','off');
    set(get(gca,'yaxis'),'MinorTick','off');
    xlabel('Frequency (Hz)');
    % ylim([-10,20]);
    % ylabel('Power (dB)');
    ylabel('Power (dB)')
    ylim([-4,10]);
    gcaformat;

    fill([8,15,15,8],[8.75,8.75,9,9],clrs(1,:),'LineStyle','none');
    % text(10.^((log10(15)+log10(8))/2),9.5,'*','FontSize',7,'HorizontalAlignment','center');
    fill([15,30,30,15],[8.75,8.75,9,9],clrs(2,:),'LineStyle','none');
    % text(10.^((log10(15)+log10(30))/2),9.5,'*','FontSize',7,'HorizontalAlignment','center');
    fill([1,4,4,1],[8.75,8.75,9,9],clrs(3,:),'LineStyle','none');
    % text(10.^((log10(0.5)+log10(4))/2),9.5,'p=0.91','FontSize',7,'HorizontalAlignment','center');

    text(10.^((log10(15)+log10(8))/2),10,'\alpha*','FontSize',7,'HorizontalAlignment','center','color',clrs(1,:));
    text(10.^((log10(15)+log10(30))/2),10,'\beta*','FontSize',7,'HorizontalAlignment','center','color',clrs(2,:));
    text(10.^((log10(1)+log10(4))/2),10,'\delta (n.s.)','FontSize',7,'HorizontalAlignment','center','color',clrs(3,:));