dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';
load(fullfile(dataFolder,'data_time_information.mat'));
t0 = timeInfo.infusion_onset-timeInfo.object_drop;


deltaIdx = find(and(freq>=0.5,freq<4));
alphaIdx = find(and(freq>=8,freq<15));
betaIdx = find(and(freq>=15,freq<30));
aligned.alpha_syn = squeeze(nanmean(aligned.syn(alphaIdx,:,:)));
aligned.beta_syn = squeeze(nanmean(aligned.syn(betaIdx,:,:)));
aligned.delta_syn = squeeze(nanmean(aligned.syn(deltaIdx,:,:)));
aligned.alpha_pre = squeeze(nanmean(aligned.pre(alphaIdx,:,:)));
aligned.beta_pre = squeeze(nanmean(aligned.pre(betaIdx,:,:)));
aligned.delta_pre = squeeze(nanmean(aligned.pre(deltaIdx,:,:)));

for i = 1:14
    % t1 = -rescaled.time*t0(i);
    % aligned.alpha_pre(:,i) = interp1(t1,rescaled.alpha_pre(:,i),aligned.time,'linear');
    % aligned.beta_pre(:,i) = interp1(t1,rescaled.beta_pre(:,i),aligned.time,'linear');
    % aligned.delta_pre(:,i) = interp1(t1,rescaled.delta_pre(:,i),aligned.time,'linear');

    % aligned.alpha_syn(:,i) = interp1(t1,rescaled.alpha_syn(:,i),aligned.time,'linear');
    % aligned.beta_syn(:,i) = interp1(t1,rescaled.beta_syn(:,i),aligned.time,'linear');
    % aligned.delta_syn(:,i) = interp1(t1,rescaled.delta_syn(:,i),aligned.time,'linear');

    aligned.delta_syn(:,i) = aligned.delta_syn(:,i)-nanmean(aligned.delta_syn(aligned.time<=t0(i),i));
    aligned.alpha_syn(:,i) = aligned.alpha_syn(:,i)-nanmean(aligned.alpha_syn(aligned.time<=t0(i),i));
    aligned.beta_syn(:,i) = aligned.beta_syn(:,i)-nanmean(aligned.beta_syn(aligned.time<=t0(i),i));


    aligned.delta_syn(:,i) = fillgaps(aligned.delta_syn(:,i),5e3);
    aligned.alpha_syn(:,i) = fillgaps(aligned.alpha_syn(:,i),1024);
    aligned.beta_syn(:,i) = fillgaps(aligned.beta_syn(:,i),512);

    aligned.delta_pre(:,i) = fillgaps(aligned.delta_pre(:,i),5e3);
    aligned.alpha_pre(:,i) = fillgaps(aligned.alpha_pre(:,i),1024);
    aligned.beta_pre(:,i) = fillgaps(aligned.beta_pre(:,i),512);
end

% Not enough data points for computing baseline
aligned.delta_syn(:,3) = nan*aligned.delta_syn(:,3);

blue = clrsPT.qualitative_CM.blue;
clrs = clrsPT.lines(3);
clrs(2:3,:) = flip(clrs(2:3,:));

% figureNB(30,10);
% for i = 1:14
%     subplot(2,7,i);
%     plot(aligned.time,aligned.delta_syn(:,i))
%     line(get(gca,'xlim'),[0,0],'color','k')
%     ylim([-10,10])
% end


aligned.alpha_pre = smoothdata(aligned.alpha_pre,'movmedian',40);
aligned.beta_pre = smoothdata(aligned.beta_pre,'movmedian',40);
aligned.delta_pre = smoothdata(aligned.delta_pre,'movmedian',40);
aligned.alpha_syn = smoothdata(aligned.alpha_syn,'movmedian',40);
aligned.beta_syn = smoothdata(aligned.beta_syn,'movmedian',40);
aligned.delta_syn = smoothdata(aligned.delta_syn,'movmedian',40);

figureNB(7,6);
subplot(3,2,1);
    plotwitherror(aligned.time,aligned.alpha_pre,'CI','LineWidth',1,'color',clrs(1,:));
    xlim([-240,60]);
    ylim([-2.75,10])
    text(-220,10,'\alpha (8-15 Hz)','FontSize',6,'color',clrs(1,:),'VerticalAlignment','top');
    ylabel('Power (dB)')
    line([-300,60],[0,0],'color','k');
subplot(3,2,3);
    plotwitherror(aligned.time,aligned.beta_pre,'CI','LineWidth',1,'color',clrs(2,:));
    xlim([-240,60]);
    ylim([-2.5,7])
    text(-220,7,'\beta (15-30 Hz)','FontSize',6,'color',clrs(2,:),'VerticalAlignment','top');
    ylabel('Power (dB)')
    line([-300,60],[0,0],'color','k');
subplot(3,2,5);
    plotwitherror(aligned.time,aligned.delta_pre,'CI','LineWidth',1,'color',clrs(3,:));
    xlim([-240,60]);
    ylim([-3.5,13])
    text(-220,13,'\delta (0.5-4 Hz)','FontSize',6,'color',clrs(3,:),'VerticalAlignment','top');
    ylabel('Power (dB)')
    line([-300,60],[0,0],'color','k');
subplot(3,2,2);
    plotwitherror(aligned.time,aligned.alpha_syn,'CI','LineWidth',1,'color',clrs(1,:));
    xlim([-240,60]);
    ylabel('Power (dB)')
    line([-300,60],[0,0],'color','k');
    ylim([-2,7]);
subplot(3,2,4);
    plotwitherror(aligned.time,aligned.beta_syn,'CI','LineWidth',1,'color',clrs(2,:));
    xlim([-240,60]);
    ylabel('Power (dB)')
    line([-300,60],[0,0],'color','k');
    % ylim([-5,10])
subplot(3,2,6);
    plotwitherror(aligned.time,aligned.delta_syn,'CI','LineWidth',1,'color',clrs(3,:));
    xlim([-240,60]);
    ylabel('Power (dB)')
    line([-300,60],[0,0],'color','k');
    ylim([-4,7])
gcaformat(gcf)

