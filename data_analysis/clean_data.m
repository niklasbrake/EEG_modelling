function clean_data2(fileName);

pathName = 'E:\Research_Projects\004_Propofol\data\experiments\rat_iso\raw';
load(fullfile(pathName,fileName));
t = [1:length(S)]'/3030;
S2 = resample(S,512,3030);
t2 = [1:length(S2)]'/512;

figureNB;
for i = 1:6
    [freq,time,psd] = eegfft(t2,S2(:,i),2,0);
    subplot(3,2,i);
    imagesc(time,freq,log10(psd));
    ylim([min(freq),50]);
    xlim([time(1),time(end)]);
    for i = 1:length(TimeInfo)
        line(TimeInfo(i)*[1,1],get(gca,'ylim'),'color','r','LineWidth',1);
    end
    set(gca,'Clim',[2,8])
    axis xy;
    gcaformat   
    drawnow;
end

fig = figureNB;
h2 = plot(freq,psd(:,1),'color',[0.6,0.6,0.6]);
hold on;
h1 = plot(freq,psd(:,1),'k');
xlim([0.5,256]);
ylim(10.^[-4,8]);
set(gca,'xscale','log');
set(gca,'yscale','log');

idcs = ones(size(psd,2),1);
for i = 1:size(psd,2)
    y = input(sprintf('time point %d',time(i)));
    if(~isempty(y))
        idcs(i) = 0;
    end
    temp = max(i-10,1):i;
    history = intersect(temp,find(idcs(temp)));
    h2.YData = nanmean(psd(:,history),2);
    h1.YData = psd(:,i);
end


data = ax.Children;
removedData = isnan(data(end).YData);
I = getIntervals(removedData);
tRemoved = t2(I);
if(length(tRemoved)==2)
    tRemoved = tRemoved(:)';
end

close(fig);

for i = 1:size(tRemoved,1)
    idcs = find(and(t>=tRemoved(i,1),t<=tRemoved(i,2)));
    S(idcs,:) = nan;
end

fig = figureNB(38,12);
fig.Position(1) = 0.6;
fig.Position(2) = 8.7;
ax = axes('Position',[0.02,0.1,0.96,0.9]);
plot(t,S+linspace(0,5e5,7));
for i = 1:length(TimeInfo)
    line(TimeInfo(i)*[1,1],get(gca,'ylim'),'color','k');
end
set(get(gca,'yaxis'),'visible','off');
xlim([t(1),t(end)]);
gcaformat;
drawnow;

[~,fileName] = fileparts(fileName);
save(fullfile(pathName,[fileName '_cleaned.mat']),'t','S','TimeInfo','Info1','Info2');
disp('saved.');
