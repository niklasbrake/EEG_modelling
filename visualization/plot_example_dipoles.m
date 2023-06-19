
load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis\sample1455\simulation\simulation_data.mat')

figureNB(4.55,4);
for i = 1:3
    plot(time,10*i+dipoles(:,i),'color','k');
    hold on;
    line([-100,-50],[10*i,10*i],'color','k','LineWidth',1)
end
xlim([-100,time(end)]);
xlim([-100,1e3])


line([400,600],[2,2],'LineWidth',1.5,'color','k')
line([375,375],[2,7],'LineWidth',1.5,'color','k')
text(355,4.5,['5 nA' char(956) 'm'],'FontSize',6,'HorizontalAlignment','right')
text(500,0,'200 ms','FontSize',6,'HorizontalAlignment','center')

text(-110,10,'P_x','FontSize',7,'HorizontalAlignment','right')
text(-110,20,'P_y','FontSize',7,'HorizontalAlignment','right')
text(-110,30,'P_z','FontSize',7,'HorizontalAlignment','right')