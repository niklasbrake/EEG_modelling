for i = 1:14
    X = TimeDomainAligned(:,2,i);
    X = fillgaps(X,1e4);


    figureNB(14,8);
    subplot(2,2,3);
        plot(Time,X,'k')
        yl = max(abs(X));
        xlim([-300,60])
        ylim([-1.1,1.1]*yl);
        xlabel('Time (s)');
        ylabel(['Voltage (' char(956) 'V)'])
        gcaformat
        hold on;
        fill([-50,50,50,-50],[-yl,-yl,yl,yl],'k','EdgeColor','r','FaceAlpha',0,'LineWidth',1)
        line([0,0],[-yl,yl],'color','r')
    subplot(2,2,1)
        plot(Time,X,'k','LineWidth',0.25);
        hold on;
        fill([-50,50,50,-50],[-yl,-yl,yl,yl],'k','EdgeColor','r','FaceAlpha',0,'LineWidth',1)
        xlim([-50,50])
        ylim([-yl,yl])
        axis off;
        line([0,0],[-yl,yl],'color','r')
        title('Broadband')



    [b0,a0] = butter(3,0.5/1024,'high');
    [b1,a1] = butter(3,4/1024,'low');
    X0 = filtfilt(b0,a0,filtfilt(b1,a1,X));

    subplot(2,2,4);
        plot(Time,X0,'k')
        yl = max(abs(X0));
        xlim([-300,60])
        ylim([-1.1,1.1]*yl);
        xlabel('Time (s)');
        ylabel(['Voltage (' char(956) 'V)'])
        gcaformat
        hold on;
        fill([-50,50,50,-50],[-yl,-yl,yl,yl],'k','EdgeColor','r','FaceAlpha',0,'LineWidth',1)
        line([0,0],[-yl,yl],'color','r')
    subplot(2,2,2)
        plot(Time,X0,'k','LineWidth',0.5);
        hold on;
        fill([-50,50,50,-50],[-yl,-yl,yl,yl],'k','EdgeColor','r','FaceAlpha',0,'LineWidth',1)
        xlim([-50,50])
        ylim([-yl,yl])
        axis off;
        line([0,0],[-yl,yl],'color','r')
        title('Bandpower (0.5-4 Hz)')
end



[b0,a0] = butter(3,8/1024,'high');
[b1,a1] = butter(3,15/1024,'low');
X0 = filtfilt(b0,a0,filtfilt(b1,a1,X));

figureNB(14,8);
subplot(2,1,2);
    plot(Time,X0,'k')
    yl = max(abs(X0));
    xlim([-300,60])
    ylim([-1.1,1.1]*yl);
    xlabel('Time (s)');
    ylabel(['Voltage (' char(956) 'V)'])
    gcaformat
    hold on;
    fill([-50,50,50,-50],[-yl,-yl,yl,yl],'k','EdgeColor','r','FaceAlpha',0,'LineWidth',1)
    line([0,0],[-yl,yl],'color','r')
subplot(2,1,1)
    plot(Time,X0,'k','LineWidth',0.1);
    hold on;
    fill([-50,50,50,-50],[-yl,-yl,yl,yl],'k','EdgeColor','r','FaceAlpha',0,'LineWidth',1)
    xlim([-50,50])
    ylim([-yl,yl])
    axis off;
    line([0,0],[-yl,yl],'color','r')
    title('Bandpower (8-15 Hz)')


[b0,a0] = butter(3,20/1024,'high');
[b1,a1] = butter(3,35/1024,'low');
X0 = filtfilt(b0,a0,filtfilt(b1,a1,X));

figureNB(14,8);
subplot(2,1,2);
    plot(Time,X0,'k')
    yl = max(abs(X0));
    xlim([-300,60])
    ylim([-1.1,1.1]*yl);
    xlabel('Time (s)');
    ylabel(['Voltage (' char(956) 'V)'])
    gcaformat
    hold on;
    fill([-50,50,50,-50],[-yl,-yl,yl,yl],'k','EdgeColor','r','FaceAlpha',0,'LineWidth',1)
    line([0,0],[-yl,yl],'color','r')
subplot(2,1,1)
    plot(Time,X0,'k','LineWidth',0.1);
    hold on;
    fill([-50,50,50,-50],[-yl,-yl,yl,yl],'k','EdgeColor','r','FaceAlpha',0,'LineWidth',1)
    xlim([-50,50])
    ylim([-yl,yl])
    axis off;
    line([0,0],[-yl,yl],'color','r')
    title('Bandpower (20-35 Hz)')