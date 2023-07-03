load(fullfile(dataFolder,'data_time_information.mat'));
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

for i = 1:14
    t1 = -rescaled.time*t0(i);
    aligned.alpha_pre(:,i) = interp1(t1,rescaled.alpha_pre(:,i),aligned.time,'linear');
    aligned.beta_pre(:,i) = interp1(t1,rescaled.beta_pre(:,i),aligned.time,'linear');
    aligned.delta_pre(:,i) = interp1(t1,rescaled.delta_pre(:,i),aligned.time,'linear');

    aligned.alpha_syn(:,i) = interp1(t1,rescaled.alpha_syn(:,i),aligned.time,'linear');
    aligned.beta_syn(:,i) = interp1(t1,rescaled.beta_syn(:,i),aligned.time,'linear');
    aligned.delta_syn(:,i) = interp1(t1,rescaled.delta_syn(:,i),aligned.time,'linear');
end


figureNB(9,4);
ax1 = axes('Position',[0.08,0.22,0.4,0.7]);
    h(1) = plotwitherror(aligned.time,smoothdata(aligned.alpha_pre,'movmedian',3),'CI','LineWidth',1,'color',clrs(1,:));
    h(2) = plotwitherror(aligned.time,smoothdata(aligned.beta_pre,'movmedian',3),'CI','LineWidth',1,'color',clrs(2,:));
    h(3) = plotwitherror(aligned.time,smoothdata(aligned.delta_pre,'movmedian',3),'CI','LineWidth',1,'color',clrs(3,:));
    xlim([-240,60]);
    ylim([-2.5,15])
    ylabel('Power (dB)')
    xlabel('Time rel. LOC (s)')
    line([-300,60],[0,0],'color','k');
    text(-220,15,'\alpha (8-15 Hz)','FontSize',6,'Color',clrs(1,:));
    text(-220,17,'\beta (15-30 Hz)','FontSize',6,'Color',clrs(2,:));
    text(-220,19,'\delta (1-4 Hz)','FontSize',6,'Color',clrs(3,:));

ax2 = axes('Position',[0.59,0.22,0.4,0.7]);
    plotwitherror(aligned.time,smoothdata(aligned.alpha_syn,'movmedian',3),'CI','LineWidth',1,'color',clrs(1,:));
    plotwitherror(aligned.time,smoothdata(aligned.beta_syn,'movmedian',3),'CI','LineWidth',1,'color',clrs(2,:));
    plotwitherror(aligned.time,smoothdata(aligned.delta_syn,'movmedian',3),'CI','LineWidth',1,'color',clrs(3,:));
    xlim([-240,60]);
    ylabel('Detrended power (dB)')
    xlabel('Time rel. LOC (s)')
    line([-300,60],[0,0],'color','k');
    ylim([-2.5,5])
    text(-220,7.5,'\alpha (8-15 Hz)','FontSize',6,'Color',clrs(1,:));
    text(-220,8.5,'\beta (15-30 Hz)','FontSize',6,'Color',clrs(2,:));
    text(-220,9.5,'\delta (1-4 Hz)','FontSize',6,'Color',clrs(3,:));



k = 4;
A = aligned.alpha_pre;
B = aligned.beta_pre;
D = aligned.delta_pre;
aligned.p = [];
M = size(A,1);
for i = 1:floor(M/k)
    aligned.p(i,1) = signtest(nanmedian(A((i-1)*k+1:i*k,:)),0,'tail','right');
    aligned.p(i,2) = signtest(nanmedian(B((i-1)*k+1:i*k,:)),0,'tail','right');
    aligned.p(i,3) = signtest(nanmedian(D((i-1)*k+1:i*k,:)),0,'tail','right');
end

A = aligned.alpha_syn;
B = aligned.beta_syn;
D = aligned.delta_syn;
aligned.p_syn = [];
aligned.time2 = [];
for i = 1:floor(M/k)
    aligned.p_syn(i,1) = signtest(nanmedian(A((i-1)*k+1:i*k,:)),0,'tail','right');
    aligned.p_syn(i,2) = signtest(nanmedian(B((i-1)*k+1:i*k,:)),0,'tail','right');
    aligned.p_syn(i,3) = signtest(nanmedian(D((i-1)*k+1:i*k,:)),0,'tail','right');
    aligned.time2(i) = aligned.time(k*(i-1)+1);
end

dt = aligned.time(2)-aligned.time(1);

idx = [2,1,3];
axes(ax1)
    for j = 1:3
        for i = 1:length(aligned.time2)
            t1 =aligned.time2(i);
            fill([t1,t1+k*dt,t1+k*dt,t1],[0,0,2,2]*(aligned.p(i,j)<0.05)+12+2*idx(j),clrs(j,:),'LineStyle','none');
            hold on
        end
    end
    xlim([-240,60]);
    % ylim([-2.5,18])
    ylim([-8,20])
    gcaformat;
    % line([0,0],[0,3],'color','k','LineWidth',1);
    % text(0.43,18.2,'*')
    % text(0.43,16.2,'*')
    % text(0.43,14.2,'*')
axes(ax2);
    for j = 1:3
        for i = 1:length(aligned.time2)
            t1 =aligned.time2(i);
            fill([t1,t1+k*dt,t1+k*dt,t1],[0,0,1,1]*(aligned.p_syn(i,j)<0.05)+6+idx(j),clrs(j,:),'LineStyle','none');
            hold on
        end
    end
    xlim([-240,60]);
    ylim([-4,10])
    gcaformat;
    % line([0,0],[0,3],'color','k','LineWidth',1);
    % text(0.43,9.1,'*')
    % text(0.43,8.1,'*')
    % text(0.43,7.1,'*')

