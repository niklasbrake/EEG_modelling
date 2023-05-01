figureNB(8.15,8.15);
    plot(f,P.m(:,1)+P.m(:,1),'LineWidth',4,'color','w');
    hold on;
    plot(f,P.m(:,1)+P.m(:,end),'LineWidth',4,'color','y');
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylim([1e-17,1e-13])
    xlim([0.1,100])
    ylabel('log power');
    yticks([]);
    xlabel('Frequency (Hz)');
    xticks(10.^[-1:2]);
    xticklabels([0.1,1,10,100]);
    % xlim([0.1,40]); set(gca,'xscale','linear');
    gcaformat_dark;
    set(get(gca,'xaxis'),'color',[0.5,0.5,0.5]);
    set(get(gca,'yaxis'),'color',[0.5,0.5,0.5])
    set(gca,'LineWidth',2)
    set(gca,'FontSize',14);


figureNB(8.15,8.15);
    plot(f,P.EI(:,1),'LineWidth',4,'color','w');
    hold on;
    plot(f,P.EI(:,end),'LineWidth',4,'color','y');
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylim([1e-17,1e-13])
    xlim([0.1,100])
    ylabel('log power');
    yticks([]);
    xlabel('Frequency (Hz)');
    xticks(10.^[-1:2]);
    xticklabels([0.1,1,10,100]);
    % xlim([0.1,40]); set(gca,'xscale','linear');
    gcaformat_dark;
    set(get(gca,'xaxis'),'color',[0.5,0.5,0.5]);
    set(get(gca,'yaxis'),'color',[0.5,0.5,0.5])
    set(gca,'LineWidth',2)
    set(gca,'FontSize',14);

figureNB(8.15,8.15);
    plot(f,P.tau(:,1),'color','w','LineWidth',4);
    hold on;
    plot(f,P.tau(:,end),'color','y','LineWidth',4);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylim([1e-17,1e-13])
    xlim([0.1,100]);
    ylabel('log power');
    yticks([]);
    xlabel('Frequency (Hz)');
    xticks(10.^[-1:2]);
    xticklabels([0.1,1,10,100]);
    % xlim([0.1,40]); set(gca,'xscale','linear');
    gcaformat_dark;
    set(get(gca,'xaxis'),'color',[0.5,0.5,0.5]);
    set(get(gca,'yaxis'),'color',[0.5,0.5,0.5])
    set(gca,'LineWidth',2)
    set(gca,'FontSize',14);

return;
figureNB;
clrs = clrsPT.lines(7); clrs = clrs([2,3,7],:);
for i = 1:3
    plot((1:length(dipoles(:,i,1)))/2e3,dipoles(:,i,1)+40*i,'color',clrs(i,:));
    hold on;
end
line([4,5],[10,10],'color','w','LineWidth',2)
line([3.8,3.8],[10,20],'color','w','LineWidth',2)
gcaformat_dark
ylim([0,140])
