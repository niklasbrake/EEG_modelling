N = 1e5;
tmax = 4e3;

[ids,ts] = simulatespikes_det(N,0,tmax*1e-3);
[ids2,ts2] = simulatespikes_det(N,0.98,tmax*1e-3);

[N1,E1] = histcounts(ts,'BinWidth',4);
t1 = 0.5*(E1(1:end-1)+E1(2:end));
N1 = N1*1e3/4/N;

[N2,E2] = histcounts(ts2,'BinWidth',4);
t2 = 0.5*(E2(1:end-1)+E2(2:end));
N2 = N2*1e3/4/N;


figureNB(3.7,2.6);
axes('Position',[0.2,0.65,0.8,0.3]);
    plot(t1,N1,'color','k');
    set(get(gca,'xaxis'),'visible','off')
    ylim([0,4]);
    gcaformat;
    text(75,4,'m=0','FontSize',6);
axes('Position',[0.2,0.15,0.8,0.3]);
    plot(t2,N2,'color','k');
    set(get(gca,'xaxis'),'visible','off')
    ylim([0,4]);
    yl = ylabel('Mean firing rate (Hz)');
    yl.Position = [-300,5.2,-1];
    gcaformat;
    line([-750,-250]+tmax,[0,0],'color','k','LineWidth',1);
    text(75,4,'m=0.98','FontSize',6);
    text(tmax-500,-1,'500 s','FontSize',6,'HorizontalAlignment','center');