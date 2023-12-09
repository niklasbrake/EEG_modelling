load(fullfile(dataFolder,'EEG_data','data_time_information.mat'));
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

rescaled = load(fullfile(dataFolder,'EEG_data','electrode2_Cz_rescaled_time.mat'));
rescaled.psd = log10(rescaled.psd);

aligned = load(fullfile(dataFolder,'EEG_data','electrode2_Cz.mat'));

load(fullfile(dataFolder,'EEG_data','Eq6_fits','electrode2_Cz.mat'),'pars');
rescaled.synPars(1:4,:,:) = permute(pars(:,1:4,:),[2,1,3]);
freq = rescaled.freq;

[full_model,synFun] = fittingmodel('eq6');
for i = 1:14
    for j = 1:size(rescaled.psd,2)
        synFit(:,j) = synFun(freq,rescaled.synPars(1:4,j,i));
    end
    rescaled.syn_detrended(:,:,i) = rescaled.psd(:,:,i)-synFit;

    rescaled.pre(:,:,i) = rescaled.psd(:,:,i)-nanmedian(rescaled.psd(:,rescaled.time<-1,i),2);
    rescaled.syn(:,:,i) = rescaled.syn_detrended(:,:,i)-nanmedian(rescaled.syn_detrended(:,rescaled.time<-1,i),2);

    synFitAligned = interp1(-t0(i)*rescaled.time',synFit',aligned.time,'cubic')';
    aligned.syn(:,:,i) = 10* (log10(aligned.psd(:,:,i)) - synFitAligned);
    P0(:,i) = nanmean(aligned.syn(:,aligned.time<t0(i),i),2);
end

alphaIdx = find(and(freq>=8,freq<15));
rescaled.alpha_pre = 10*squeeze(nanmean(rescaled.pre(alphaIdx,:,:)));
rescaled.alpha_syn = 10*squeeze(nanmean(rescaled.syn(alphaIdx,:,:)));
aligned.alpha_pre = squeeze(nanmean(aligned.pre(alphaIdx,:,:)));
aligned.alpha_syn = squeeze(nanmean(aligned.syn(alphaIdx,:,:)));

betaIdx = find(and(freq>=15,freq<30));
rescaled.beta_pre = 10*squeeze(nanmean(rescaled.pre(betaIdx,:,:)));
rescaled.beta_syn = 10*squeeze(nanmean(rescaled.syn(betaIdx,:,:)));
aligned.beta_pre = squeeze(nanmean(aligned.pre(betaIdx,:,:)));
aligned.beta_syn = squeeze(nanmean(aligned.syn(betaIdx,:,:)));

deltaIdx = find(and(freq>=1,freq<=4));
rescaled.delta_pre = 10*squeeze(nanmean(rescaled.pre(deltaIdx,:,:)));
rescaled.delta_syn = 10*squeeze(nanmean(rescaled.syn(deltaIdx,:,:)));
aligned.delta_syn = squeeze(nanmean(aligned.syn(deltaIdx,:,:)));
aligned.delta_pre = squeeze(nanmean(aligned.pre(deltaIdx,:,:)));

ptIdx = 1;
preExample = 10.^nanmedian(rescaled.psd(:,rescaled.time<-1,ptIdx),2);
idcs = find(and(aligned.time>0,aligned.time<10));
postExample = nanmedian(aligned.psd(:,idcs,ptIdx),2);
synPre = synDetrend(freq(freq<100),preExample(freq<100),3,'eq6');
synPost = synDetrend(freq(freq<100),postExample(freq<100),3,'eq6');

% Write to Source Data files
% filename = 'E:\Research_Projects\004_Propofol\manuscript\Nature Communications\_final_submission\source_data.xlsx';
% x1 = round(log10(preExample),3,'significant');
% x2 = round(log10(postExample),3,'significant');
% T = table(freq,x1(:),x2(:));
% T.Properties.VariableNames = {'Time (s)','Baseline','post-LOC'};
% writetable(T,filename,'Sheet','Figure 9a, e','Range','B2')

% x = nanmean(aligned.pre,3);
% T = table(freq(freq<=40));
% T.Properties.VariableNames{1} = 'Frequency (Hz)';
% for i = 1:size(x,2)
%     T{:,i+1} = round(x(freq<=40,i),3,'significant');
%     T.Properties.VariableNames{i+1} = sprintf('%.1f s',aligned.time(i));
% end
% writetable(T,filename,'Sheet','Figure 9b','Range','B2')

% x = nanmean(aligned.syn-permute(P0,[1,3,2]),3);
% T = table(freq(freq<=40));
% T.Properties.VariableNames{1} = 'Frequency (Hz)';
% for i = 1:size(x,2)
%     T{:,i+1} = round(x(freq<=40,i),3,'significant');
%     T.Properties.VariableNames{i+1} = sprintf('%.1f s',aligned.time(i));
% end
% writetable(T,filename,'Sheet','Figure 9f','Range','B2')

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
    ylim([1e-4,1e4]);
    yticks(10.^[-4:4:4]);
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
    ylim([-15,20]);
    ylabel('Power (dB)')
axes('Position',[0.34,0.62,0.19,0.3]);
    imagesc(aligned.time,freq,nanmean(aligned.pre,3));
    ylim([0.5,50])
    CB = colorbar('location','eastoutside');
    CB.Label.String = 'Power (dB)';
    CB.Position(1) = 0.54;
    axis xy
    xlim([-180,60]);
    xticks([-180:60:60]);
    ylabel('Frequency (Hz)')
    xlabel('LOC-aligned time (s)')
    set(gca,'CLim',[-5,10])
ax1 = axes('Position',[0.82,0.62,0.15,0.3]);
    h(1) = plotwitherror(rescaled.time,smoothdata(rescaled.alpha_pre,'movmedian',3),'CI','LineWidth',1,'color',clrs(1,:));
    h(2) = plotwitherror(rescaled.time,smoothdata(rescaled.beta_pre,'movmedian',3),'CI','LineWidth',1,'color',clrs(2,:));
    h(3) = plotwitherror(rescaled.time,smoothdata(rescaled.delta_pre,'movmedian',3),'CI','LineWidth',1,'color',clrs(3,:));
    xlim([-1.25,0.25]);
    ylim([-2.5,15])
    ylabel('Power (dB)')
    xlabel('Rescaled time')
    xticks([-1,0]);
    xticklabels({'Infusion','LOC'})
    line([-1.5,0.5],[0,0],'color','k');
    text(-1.3,9.5,'\alpha (8-15 Hz)','FontSize',6,'Color',clrs(1,:));
    text(-1.3,7,'\beta (15-30 Hz)','FontSize',6,'Color',clrs(2,:));
    text(-1.3,12,'\delta (1-4 Hz)','FontSize',6,'Color',clrs(3,:));

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
    ylim([1e-4,1e4]);
    yticks(10.^[-4:4:4]);
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
    ylim([-15,20]);
    ylabel('Power (dB)')
axes('Position',[0.34,0.15,0.19,0.3]);
    imagesc(aligned.time,aligned.freq,nanmean(aligned.syn-permute(P0,[1,3,2]),3))
    ylim([0.5,50])
    CB = colorbar('location','eastoutside');
    CB.Label.String = 'Power (dB)';
    CB.Position(1) = 0.54;
    axis xy
    xlim([-180,60]);
    xticks([-180:60:60]);
    ylabel('Frequency (Hz)')
    xlabel('LOC-aligned time (s)')
    set(gca,'CLim',[-3,6])
ax2 = axes('Position',[0.82,0.15,0.15,0.3]);
    plotwitherror(rescaled.time,smoothdata(rescaled.alpha_syn,'movmedian',3),'CI','LineWidth',1,'color',clrs(1,:));
    plotwitherror(rescaled.time,smoothdata(rescaled.beta_syn,'movmedian',3),'CI','LineWidth',1,'color',clrs(2,:));
    plotwitherror(rescaled.time,smoothdata(rescaled.delta_syn,'movmedian',3),'CI','LineWidth',1,'color',clrs(3,:));
    xlim([-1.25,0.25]);
    ylabel('Power (dB)')
    xlabel('Rescaled time')
    xticks([-1,0]);
    xticklabels({'Infusion','LOC'})
    line([-1.5,0.5],[0,0],'color','k');
    ylim([-2.5,5])
    text(-1.3,9.5/2,'\alpha (8-15 Hz)','FontSize',6,'Color',clrs(1,:));
    text(-1.3,7/2,'\beta (15-30 Hz)','FontSize',6,'Color',clrs(2,:));
    text(-1.3,12/2,'\delta (1-4 Hz)','FontSize',6,'Color',clrs(3,:));

gcaformat(fig);

idcs = [interp1(linspace(0,1,200),1:200,linspace(0,1,50),'nearest','extrap'),interp1(linspace(0,1,200),201:400,linspace(0,1,100),'nearest','extrap')];
CM = clrsPT.iridescent(400);
colormap(flip(CM(idcs,:)));

% Write to Source Data files
% time = linspace(-0.5,1.5,200);
% T = table(time(:));
% T.Properties.VariableNames{1} = 'Time (rescaled)';
% for i = 1:14
%     T{:,i+1} = round(rescaled.delta_pre(:,i),3,'significant');
%     T.Properties.VariableNames{i+1} = sprintf('pt. %d, delta',i);
% end
% for i = 1:14
%     T{:,i+15} = round(rescaled.alpha_pre(:,i),3,'significant');
%     T.Properties.VariableNames{i+15} = sprintf('pt. %d, alpha',i);
% end
% for i = 1:14
%     T{:,i+29} = round(rescaled.beta_pre(:,i),3,'significant');
%     T.Properties.VariableNames{i+29} = sprintf('pt. %d, beta',i);
% end
% writetable(T,filename,'Sheet','Figure 9d','Range','B2')

% time = linspace(-0.5,1.5,200);
% T = table(time(:));
% T.Properties.VariableNames{1} = 'Time (rescaled)';
% for i = 1:14
%     T{:,i+1} = round(rescaled.delta_syn(:,i),3,'significant');
%     T.Properties.VariableNames{i+1} = sprintf('pt. %d, delta',i);
% end
% for i = 1:14
%     T{:,i+15} = round(rescaled.alpha_syn(:,i),3,'significant');
%     T.Properties.VariableNames{i+15} = sprintf('pt. %d, alpha',i);
% end
% for i = 1:14
%     T{:,i+29} = round(rescaled.beta_syn(:,i),3,'significant');
%     T.Properties.VariableNames{i+29} = sprintf('pt. %d, beta',i);
% end
% writetable(T,filename,'Sheet','Figure 9h','Range','B2')


k = 5;
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

alpha = 0.05;

idx = [2,1,3];
axes(ax1)
    for j = 1:3
        for i = 1:length(rescaled.time2)
            t1 =rescaled.time2(i);
            fill([t1,t1+k*dt,t1+k*dt,t1],[0,0,2,2]*(rescaled.p(i,j)<alpha)+12+2*idx(j),clrs(j,:),'LineStyle','none');
            hold on
        end
    end
    xlim([-1.4,0.4]);
    ylim([-8,20])
    xticks([-1,-0.5,0]);
    xticklabels({'Infusion','','LOC'});
    text(0.43,18.2,'*')
    text(0.43,16.2,'*')
    text(0.43,14.2,'*')
axes(ax2);
    for j = 1:3
        for i = 1:length(rescaled.time2)
            t1 =rescaled.time2(i);
            fill([t1,t1+k*dt,t1+k*dt,t1],[0,0,1,1]*(rescaled.p_syn(i,j)<alpha)+6+idx(j),clrs(j,:),'LineStyle','none');
            hold on
        end
    end
    xlim([-1.4,0.4]);
    ylim([-4,10])
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
A = nanmean(post(alphaIdx,:)); fprintf('\\alpha: p=%f\n',signtest(A,0,'tail','right'));
B = nanmean(post(betaIdx,:)); fprintf('\\beta: p=%f\n',signtest(B,0,'tail','right'));
D = nanmean(post(deltaIdx,:)); fprintf('\\delta: p=%f\n',signtest(D,0,'tail','right'));
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
    ylim([-8,20]);
    gcaformat;

    fill([8,15,15,8],[18.5,18.5,19,19]-1,clrs(1,:),'LineStyle','none');
    fill([15,30,30,15],[18.5,18.5,19,19]-1,clrs(2,:),'LineStyle','none');
    fill([1,4,4,1],[18.5,18.5,19,19]-1,clrs(3,:),'LineStyle','none');

    text(10.^((log10(15)+log10(8))/2),20,'\alpha*','FontSize',7,'HorizontalAlignment','center','color',clrs(1,:));
    text(10.^((log10(15)+log10(30))/2),20,'\beta*','FontSize',7,'HorizontalAlignment','center','color',clrs(2,:));
    text(10.^((log10(1)+log10(4))/2),20,'\delta*','FontSize',7,'HorizontalAlignment','center','color',clrs(3,:));

for i = 1:14
    tt = -t0(1)*rescaled.time;
    idcs = find(and(tt>=-10,tt<=0));
    post_syn(:,i) = 10*nanmedian(rescaled.syn(:,idcs,i),2);
end
A_syn = nanmean(post_syn(alphaIdx,:)); fprintf('\\alpha: p=%f\n',signtest(A_syn,0,'tail','right'));
B_syn = nanmean(post_syn(betaIdx,:)); fprintf('\\beta: p=%f\n',signtest(B_syn,0,'tail','right'));
D_syn = nanmean(post_syn(deltaIdx,:)); fprintf('\\delta: p=%f\n',signtest(D_syn,0,'tail','right'));
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
    ylabel('Power (dB)')
    ylim([-4,10]);
    gcaformat;

    fill([8,15,15,8],[8.75,8.75,9,9],clrs(1,:),'LineStyle','none');
    fill([15,30,30,15],[8.75,8.75,9,9],clrs(2,:),'LineStyle','none');
    fill([1,4,4,1],[8.75,8.75,9,9],clrs(3,:),'LineStyle','none');

    text(10.^((log10(15)+log10(8))/2),10,'\alpha*','FontSize',7,'HorizontalAlignment','center','color',clrs(1,:));
    text(10.^((log10(15)+log10(30))/2),10,'\beta*','FontSize',7,'HorizontalAlignment','center','color',clrs(2,:));
    text(10.^((log10(1)+log10(4))/2),10,'\delta (n.s.)','FontSize',7,'HorizontalAlignment','center','color',clrs(3,:));

% Write to Source Data files
% T = table(rescaled.freq);
% T.Properties.VariableNames{1} = 'Frequency (Hz)';
% for i = 1:14
%     T{:,i+1} = round((post(:,i)),3,'significant');
%     T.Properties.VariableNames{i+1} = sprintf('pt. %d',i);;
% end
% writetable(T,filename,'Sheet','Figure 9c','Range','B2')

% T = table(rescaled.freq);
% T.Properties.VariableNames{1} = 'Frequency (Hz)';
% for i = 1:14
%     T{:,i+1} = round((post_syn(:,i)),3,'significant');
%     T.Properties.VariableNames{i+1} = sprintf('pt. %d',i);;
% end
% writetable(T,filename,'Sheet','Figure 9g','Range','B2')
