
% load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis\sample1455\simulation\simulation_data.mat')
load(fullfile(dataFolder,'example_simulation','simulation','simulation_data.mat'))
load(fullfile(dataFolder,'example_simulation','model.mat'))

% Resample to 1KHz
dp = resample(dipoles,1,16);
time = time(1:16:end);

figureNB(4.55,4);
for i = 1:3
    plot(time,20*i+dp(:,i),'color','k');
    hold on;
    line([-100,-50],[20*i,20*i],'color','k','LineWidth',1)
end
xlim([-100,1e3])
axis off


line([400,600],[0,0],'LineWidth',1.5,'color','k')
line([375,375],[0,10],'LineWidth',1.5,'color','k')
text(355,4.5,['10 nA' char(956) 'm'],'FontSize',6,'HorizontalAlignment','right')
text(500,-4,'200 ms','FontSize',6,'HorizontalAlignment','center')

text(-110,20,'P_x','FontSize',7,'HorizontalAlignment','right')
text(-110,40,'P_y','FontSize',7,'HorizontalAlignment','right')
text(-110,60,'P_z','FontSize',7,'HorizontalAlignment','right')



[ids,ts,ei] = network.getprenetwork(network.spikingFile);
nI = sum(ei);
n = length(ei);
figureNB(3.1,2.8);
    R = raster(ids(ei(ids)==0),ts(ei(ids)==0),gcf);
    R.Color = 'r'; R.LineWidth = 0.5;
    hold on
    R = raster(ids(ei(ids)==1),ts(ei(ids)==1),gcf);
    R.Color = 'b'; R.LineWidth = 0.5;
    ylim([1,length(ei)])
    yticks([1,n-nI,n])
    xticks([1,2e3]);
    yticklabels({1,'',n});
    xlabel('Time (ms)');
    yl = ylabel('Synapse IDs');
    yl.Position = [-450,5350,-1];
    set(get(gca,'yaxis'),'visible','on')
    gcaformat;