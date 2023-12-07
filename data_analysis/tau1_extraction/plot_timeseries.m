eeg_example = load(fullfile(dataFolder,'data_sample_time_series.mat'));
psd = [];

load(fullfile(dataFolder,'data_time_information.mat'));
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

load(fullfile(dataFolder,'data_Cz_multitaper_meanRef.mat'));

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

CM = clrsPT.iridescent(1e3);
black = CM(end,:);

marsh_model_analysis;

fig = figureNB(8.9,10);
axes('Position',[0.067,0.76,0.91,0.2]);
    plot(downsample(eeg_example.time,50),downsample(eeg_example.timedomain,50),'LineWidth',0.2,'color',black); hold on;
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
axes('Position',[0.067,0.69,0.91,0.05]);
    x = t{1}; y = p{1};
    idcs = find(y<0.01);
    y(idcs)=nan;
    plot(x,y,'color',red,'LineWidth',1);
    xlim([eeg_example.time(1)-10,eeg_example.time(end)]);
    yl = get(gca,'ylim');
    axis off;

x1 = eeg_example.time(1)-10;
x2 = eeg_example.time(end);
a = (x(idcs(end))-x1)/(x2-x1)
axes('Position',[0.067+0.91*a-0.02,0.69,0.01,0.05]);
    set(get(gca,'xaxis'),'visible','off')
    set(gca,'color','none');
    gcaformat
    ylim(yl);
    yticks([0,1,2])
    gcaformat