
dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';

load(fullfile(dataFolder,'data_time_information.mat'));
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

psd = [];
aligned = load(fullfile(dataFolder,'data_aligned_detrended.mat'));
rescaled = load(fullfile(dataFolder,'data_rescaled_detrended.mat'));
freq = rescaled.freq;

figureNB(21,4.3);
for i = 1:14
    subplot(2,14,i);
        y = 10.^nanmedian(rescaled.psd(:,rescaled.time<-1,i),2);
        % y0 = 10.^synFun(freq,p0(:,i));
        pp = p0(1:4,i); pp(2) = min(pp(2),50e-3);
        [px,synFun,full_model] = synDetrend(freq(freq<100),y(freq<100),1,'exp2',[1e-2,4e-3,-10,4]);
        y0 = 10.^synFun(freq,px);
        plot(freq,y,'k','LineWidth',1); hold on;
        h0 = plot(freq,y0,'LineWidth',1,'color','b');  hold on
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        xlim([0.5,100]);
        xticks([0.5,5,50]);
        xticklabels([0.5,5,50]);
        set(get(gca,'xaxis'),'MinorTick','off');
        set(get(gca,'yaxis'),'MinorTick','off');
        xlabel('Frequency (Hz)');
        ylim([1e-2,1e4]);
        yticks(10.^[-2:2:4]);
        ylabel(['PSD (' char(956) 'V^2/Hz)']);
        xlabel('');
        xticklabels({});
        text(1,5e3,['Pt. ' int2str(i)],'FontSize',6)
        text(0.7,5e-2,[num2str(p0(1,i)*1e3,3) ' ms'],'FontSize',6,'color','b')
        gcaformat
    subplot(2,14,i+14);
        h1 = plot(freq,10*log10(y./y0),'color','k','LineWidth',1);
        set(gca,'xscale','log');
        xlim([0.5,100]);
        xticks([0.5,5,50]);
        xticklabels([0.5,5,50]);
        ylim([-5,20])
        line([0.5,100],[0,0],'color','r','linestyle','--');
        set(get(gca,'xaxis'),'MinorTick','off');
        set(get(gca,'yaxis'),'MinorTick','off');
        % xlabel('Frequency (Hz)');
        ylabel('Power (dB)');
        gcaformat
    if(i~=1 && i~=15)
        % set(get(gca,'yaxis'),'visible','off');
        subplot(2,14,i);
        ylabel('');
        yticklabels({});
        subplot(2,14,i+14);
        ylabel('');
        yticklabels({});
    end
    drawnow;
    continue;;
    fprintf('p = [%.3f,%.3f,%.3f,%.3f]\n',p0(1:4,i));
    while(true)
        x = input('p1 = ');
        if(~isnumeric(x))
            break;
        else
            p0(1:4,i) = x(1:4);
            h0.YData = 10.^synFun(freq,p0(:,i));
            h1.YData = 10*log10(y./10.^synFun(freq,p0(:,i)));
        end
    end
end

figureNB(21,4.3);
% figureNB(40,8);
for i = 1:14
    subplot(2,14,i);
        idx = find(and(rescaled.time*-t0(i)>0-10,rescaled.time*-t0(i)<=0));
        y = 10.^nanmedian(rescaled.psd(:,idx,i),2);
        % y0 = 10.^synFun(freq,p1(:,i));
        pp = p1(1:4,i); pp(2) = min(pp(2),50e-3);
        [px,synFun,full_model] = synDetrend(freq(freq<100),y(freq<100),2,'exp2',pp);
        y0 = 10.^synFun(freq,px);
        plot(freq,y,'k','LineWidth',1); hold on;
        h0 = plot(freq,y0,'LineWidth',1,'color','b');  hold on
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        xlim([0.5,100]);
        xticks([0.5,5,50]);
        xticklabels([0.5,5,50]);
        set(get(gca,'xaxis'),'MinorTick','off');
        set(get(gca,'yaxis'),'MinorTick','off');
        xlabel('Frequency (Hz)');
        ylim([1e-4,3e4]);
        yticks(10.^[-4:4:4]);
        ylabel(['PSD (' char(956) 'V^2/Hz)']);
        xlabel('');
        xticklabels({});
        text(1,10e3,['Pt. ' int2str(i)],'FontSize',6)
        text(0.7,1e-3,[num2str(px(1)*1e3,3) ' ms'],'FontSize',6,'color','b')
        gcaformat
    subplot(2,14,i+14);
        h1 = plot(freq,10*log10(y./y0),'color','k','LineWidth',1);
        set(gca,'xscale','log');
        xlim([0.5,100]);
        xticks([0.5,5,50]);
        xticklabels([0.5,5,50]);
        ylim([-5,20]);
        line([0.5,100],[0,0],'color','r','linestyle','--');
        set(get(gca,'xaxis'),'MinorTick','off');
        set(get(gca,'yaxis'),'MinorTick','off');
        % xlabel('Frequency (Hz)');
        ylabel('Power (dB)');
        gcaformat
    if(i~=1 && i~=15)
        % set(get(gca,'yaxis'),'visible','off');
        subplot(2,14,i);
        ylabel('');
        yticklabels({});
        subplot(2,14,i+14);
        ylabel('');
        yticklabels({});
    end
    drawnow;
    continue;;
    fprintf('p = [%.3f,%.3f,%.3f,%.3f]\n',p1(1:4,i));
    while(true)
        x = input('p1 = ');
        if(~isnumeric(x))
            break;
        else
            p1(1:4,i) = x(1:4);
            h0.YData = 10.^synFun(freq,p1(:,i));
            h1.YData = 10*log10(y./10.^synFun(freq,p1(:,i)));
        end
    end
end


load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\raw\time_series_all_channels.mat')

pt = 13;
t0 = -113;
idcs = find(and(Time>t0,Time<=t0+10));
t = Time(idcs);
Y = TimeDomainAligned(idcs,2,pt);

[b,a] = butter(7,55/1024*2,'low');

[freq,time,psd] = eegfft(t,Y,2,0);

psd(and(freq>55,freq<65),:) = nan;
for i = 1:size(psd,2)
    psd(:,i) = fillgaps(psd(:,i),5);
end



figureNB(12,6);
subplot(2,1,1);
    plot(t,filtfilt(b,a,Y),'color','k')
    hold on;
    % plot(t,filtfilt(b,a,Y),'r','LineWidth',2)
    ylim([-50,60]);
    xlim([t0,t0+10]);
    xlabel('LOC-aligned time (s)')
    ylabel(['EEG (' char(956) 'V)'])
    set(gca,'FontSize',7)
    gcaformat;
    box on;

for i = 1:5
    t1 = t0+(i-1)*2;
    fill([t1,t1+2,t1+2,t1],[-50,-50,60,60],'r','FaceAlpha',0.1,'EdgeColor','r','LineWidth',1);
end

for i = 1:5
    ax = subplot(2,5,5+i);
    plot(freq,psd(:,i),'LineWidth',1,'color','k')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.5,100]);
    ylim([1e-2,1e2]);
    P(:,i) = synDetrend(freq(freq<100),psd(freq<100,i),2,'exp2',[0.015,4e-3,-10,4]);
        hold on;
    plot(freq,10.^synFun(freq,P(:,i)),'LineWidth',1,'color','r')
    drawnow
    % xax = get(gca,'xaxis');
    % xax.Color = 'r';
    % yax = get(gca,'yaxis');
    % yax.Color = 'r';
    xlabel('Frequency (Hz)');
    set(gca,'FontSize',7)
    gcaformat;
    box on;
    ax.Position(2) = 0.127;
    ax.Position(4) = 0.324;
    if(i>1)
        yticklabels({})
    else
        ylabel(['PSD (' char(956) 'V^2/Hz)']);
    end
end

tcov = linspace(0,0.300,1e3);
ecov = exp(-tcov/mean(P(1,:)))-exp(-tcov/mean(P(2,:)));
ecov = ecov/sum(ecov);

dt = 1e-3;
tcov = 0:dt:0.2;
ecov = exp(-tcov/mean(P(1,:)))-exp(-tcov/mean(P(2,:)));
ecov = ecov/sum(ecov);
N = randn(10/dt,1);


t0 = dt*(1:length(N))';
Yp = filter(b,a,60*(filter(ecov,1,N)+0.1*randn(size(N))+0.2*sin(2*pi*t0(:)*10)));
subplot(3,1,3);
    plot(t0,Yp,'color','k');
    ylim([-50,60]);