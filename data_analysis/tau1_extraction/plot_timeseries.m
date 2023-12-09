load(fullfile(dataFolder,'EEG_data','data_sample_time_series.mat'),'time','eeg','ptID');

load(fullfile(dataFolder,'EEG_data','data_time_information.mat'));
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

red = clrsPT.qualitative_CM.red;

CM = clrsPT.iridescent(1e3);
black = CM(end,:);

marsh_model_analysis;

fig = figureNB(8.9,10);
axes('Position',[0.067,0.76,0.91,0.2]);
    plot(time,eeg,'LineWidth',0.2,'color',black); hold on;
    xlim([time(1)-10,time(end)]);
    ylim([-100,100]);
    line([time(1)-10,time(1)-10],[-25,25],'color','k','linewidth',1);
    line([time(1),time(1)+30],[-50,-50],'color','k','linewidth',1);
    text(time(1)-15,0,['50 ' char(956) 'V'],'VerticalAlignment','bottom','HorizontalAlignment','center','fontsize',7,'color','k','Rotation',90);
    text(time(1)+15,-55,'30 s','VerticalAlignment','top','HorizontalAlignment','center','fontsize',7,'color','k');
    scatter(0,90,10,'vk','filled');
    text(0,105,sprintf('LOC'),'FontSize',7,'VerticalAlignment','bottom','HorizontalAlignment','center','color','k');
    line([t0(ptID),0],[75,75],'color',red,'linewidth',2);
    text(-50,80,'Propofol infusion','color',red,'FontSize',7,'VerticalAlignment','bottom','HorizontalAlignment','center');
    axis off;
axes('Position',[0.067,0.69,0.91,0.05]);
    x = t{ptID}; y = p{ptID};
    idcs = find(y<0.01);
    y(idcs)=nan;
    plot(x,y,'color',red,'LineWidth',1);
    xlim([time(1)-10,time(end)]);
    yl = get(gca,'ylim');
    axis off;

x1 = time(1)-10;
x2 = time(end);
a = (x(idcs(end))-x1)/(x2-x1)
axes('Position',[0.067+0.91*a-0.02,0.69,0.01,0.05]);
    set(get(gca,'xaxis'),'visible','off')
    set(gca,'color','none');
    gcaformat
    ylim(yl);
    yticks([0,1,2])
    gcaformat