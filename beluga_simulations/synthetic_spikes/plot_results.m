% Import all data
folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\synthetic_spikes\mCombo1';
load(fullfile(folder,'model.mat'));
network = network.rebaseNetwork(folder);
load(fullfile(folder,'simulation','simulation_data.mat'));
N = network.getsynapsecount;

% Get synapse spikes times
cons = csvread(fullfile(network.postNetwork,'connections.csv'));
[ids,ts,ei] = network.getprenetwork(network.spikingFile);
[~,I] = sort(cons(:,3));
cons = cons(I,:);

% Get synapse locations
mTypeSegmentationData = fullfile(network_simulation_beluga.resourceFolder,'cortical_column_Hagen','morphology_segmentations.mat');
load(mTypeSegmentationData)
for i = 1:length(network.neurons)
    mData{i} = nrnSegs.(network.neurons(i).mType);
    mData{i}.pos = [mean(mData{i}.x,2), ...
    mean(mData{i}.y,2), ...
    mean(mData{i}.z,2)];
end
P = arrayfun(@(i,j) mData{i}.pos(j,:),cons(:,1),cons(:,2),'UniformOutput',false);
P = cat(1,P{:});
[thet,phi] = cart2sph(P(:,1),P(:,2),P(:,3));

% Order syanpses by space
[~,J] = sort(thet);
thet = thet(J);
phi = phi(J);
synSphCoord = [phi(:),thet(:)];
[~,I] = sort(J);
ids = I(ids);

% nrnID = cons(I,1);
% [nrnID,J] = sort(nrnID);
% [~,I] = sort(J);
% ids = J(ids);

[ids_spaceshuffle,ts_spaceshuffle] = network.getprenetwork(fullfile(network.preNetwork,'spikeTimes_spaceshuffle.csv'));
[ids_timeshuffle,ts_timeshuffle] = network.getprenetwork(fullfile(network.preNetwork,'spikeTimes_timeshuffle.csv'));

ts0 = cell(N,1);
ts0_timeshuffle = cell(N,1);
ts0_spaceshuffle = cell(N,1);
for i = 1:N
    ts0{i} = ts(ids==i);
    ts0_timeshuffle{i} = ts_timeshuffle(ids_timeshuffle==i);
    ts0_spaceshuffle{i} = ts_spaceshuffle(ids_spaceshuffle==i);
end

idcs = randi(N,1e6,2);
C = zeros(length(idcs),1);
C_timeshuffle= zeros(length(idcs),1);
C_spaceshuffle = zeros(length(idcs),1);
havDist = zeros(length(idcs),1);
hav = @(x) (1-cos(x))/2;
hav_d = @(p1,x) hav(p1(1)-x(:,1))+(1-hav(x(:,1)-p1(1))-hav(p1(1)+x(:,1))).*hav(p1(2)-x(:,2));
for i = 1:length(idcs)
    i0 = idcs(i,1);
    j0 = idcs(i,2);
    havDist(i) = hav_d(synSphCoord(i0,:),synSphCoord(j0,:));
    ti = ts0{i0};
    tj = ts0{j0};
    dif = abs(ti-tj');
    if(length(dif(:))==0)
        C(i) = 0;
    else
        C(i) = 100*sum(dif(:)<2)/length(dif(:));
    end
    ti = ts0_timeshuffle{i0};
    tj = ts0_timeshuffle{j0};
    dif = abs(ti-tj');
    if(length(dif(:))==0)
        C_timeshuffle(i) = 0;
    else
        C_timeshuffle(i) = 100*sum(dif(:)<2)/length(dif(:));
    end
    ti = ts0_spaceshuffle{i0};
    tj = ts0_spaceshuffle{j0};
    dif = abs(ti-tj');
    if(length(dif(:))==0)
        C_spaceshuffle(i) = 0;
    else
        C_spaceshuffle(i) = 100*sum(dif(:)<2)/length(dif(:));
    end
end

hBins = 10.^linspace(-2,0,10);
dh = min(diff(hBins))
hEdges = [hBins(1)-dh,hBins+dh/2];

M = zeros(size(hBins));
M_timeshuffle = zeros(size(hBins));
M_spaceshuffle = zeros(size(hBins));
SD = zeros(size(hBins));
SD_timeshuffle = zeros(size(hBins));
SD_spaceshuffle = zeros(size(hBins));
for i = 1:length(hEdges)-1
    idcs = find(and(havDist>hEdges(i),havDist<hEdges(i+1)));
    idcs = randsample(idcs,1000);
    M(i) = mean(C(idcs));
    SD(i) = stderror(C(idcs));
    M_spaceshuffle(i) = mean(C_spaceshuffle(idcs));
    SD_spaceshuffle(i) = stderror(C_spaceshuffle(idcs));
    M_timeshuffle(i) = mean(C_timeshuffle(idcs));
    SD_timeshuffle(i) = stderror(C_timeshuffle(idcs));
end


% clrs = clrsPT.lines(3); clrs(2:end,:);
clrs(1,:) =  [0,0,0];
clrs(2,:) = clrsPT.qualitative_CM.cyan;
clrs(3,:) = clrsPT.qualitative_CM.teal;
figureNB(9.2,4)
subplot(1,3,1);
    fill([hBins,fliplr(hBins)],[M-1.96*SD,fliplr(M+1.96*SD)],clrs(1,:),'LineStyle','none','FaceAlpha',0.25);
    hold on;
    plot(hBins,M,'-k','LineWidth',1);
    fill([hBins,fliplr(hBins)],[M_timeshuffle-1.96*SD_timeshuffle,fliplr(M_timeshuffle+1.96*SD_timeshuffle)],clrs(2,:),'LineStyle','none','FaceAlpha',0.25);
    hold on;
    plot(hBins,M_timeshuffle,'color',clrs(2,:),'LineWidth',1)
    fill([hBins,fliplr(hBins)],[M_spaceshuffle-1.96*SD_spaceshuffle,fliplr(M_spaceshuffle+1.96*SD_spaceshuffle)],clrs(2,:),'LineStyle','none','FaceAlpha',0.25);
    hold on;
    plot(hBins,M_spaceshuffle,'color',clrs(3,:),'LineWidth',1)
    ylabel('Synchronous inputs (%)')
    xlabel('Angular seperation')
    gcaformat
    xlim([-0.05,1]);


% Get pdf of input synchronization
h = histcounts(ts,'BinEdges',0.5:network.tmax+0.5);
[P,f] = pspectrum(detrend(h,'constant'),1e3);
[synchCount,bins] = histcounts(h,'BinEdges',[-0.5:40],'Normalization','pdf');
h_timeshuffle = histcounts(ts_timeshuffle,'BinEdges',0.5:network.tmax+0.5);
[P_timeshuffle,f] = pspectrum(detrend(h_timeshuffle,'constant'),1e3);
[synchCount_timeshuffle,bins] = histcounts(h_timeshuffle,'BinEdges',[-0.5:40],'Normalization','pdf');
h_spaceshuffle = histcounts(ts_spaceshuffle,'BinEdges',0.5:network.tmax+0.5);
[P_spaceshuffle,f] = pspectrum(detrend(h_spaceshuffle,'constant'),1e3);
[synchCount_spaceshuffle,bins] = histcounts(h_spaceshuffle,'BinEdges',[-0.5:40],'Normalization','pdf');

subplot(1,3,2);
    plot(bins(1:end-1)+0.5,synchCount,'color',clrs(1,:),'LineWidth',1);
    hold on;
    plot(bins(1:end-1)+0.5,synchCount_timeshuffle,'-','color',clrs(2,:),'LineWidth',1);
    plot(bins(1:end-1)+0.5,synchCount_spaceshuffle,'--','color',clrs(3,:),'LineWidth',1);
    xlim([-2,40]);
    % set(gca,'yscale','log')
    xlabel('Synchronous spikes / ms')
    ylabel('Frequency');
    gcaformat;
subplot(1,3,3);
    plot(f,P,'color',clrs(1,:)); hold on;
    plot(f,P_timeshuffle,'-','color',clrs(2,:));
    plot(f,P_spaceshuffle,'--','color',clrs(3,:));
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xlabel('Frequency (Hz)');
    xticks([1,100])
    xticklabels([1,100])
    ylabel('PSD (spike count^2/Hz)');
    gcaformat;

t = 1:network.tmax;
figureNB(5,3.8);
axes('Position',[0.1,0.77,0.84,0.2]);
    plot(t,h/length(ei)*1e3,'k');
    xlim([0,5e3]);
    axis off;
    yl = get(gca,'ylim');
axes('Position',[0.085,0.77,0.015,0.2]);
    ylim(yl);
    set(get(gca,'xaxis'),'visible','off');
    gcaformat;
axes('Position',[0.1,0.16,0.84,0.56]);
    idcs = find(nrnID(ids)==1);
    raster(ids(idcs),ts(idcs),gcf);
    hold on;
    idcs = find(nrnID(ids)==2);
    R = raster(ids(idcs),ts(idcs),gcf);
    R.Color='r'
    xlim([0,5e3]);
    ylim([1,max(ids)]);
    yticks([1,min(ids(idcs)),max(ids)]);
    line(get(gca,'xlim'),[1,1]*min(ids(idcs)),'color','k');
    xlabel('Time (ms)');
    yticks([])
    ylabel('Synaptic input raster')
    gcaformat;
    set(get(gca,'yaxis'),'visible','on')
    box on;

folders = {fullfile(folder,'simulation','simulation_data.mat'),fullfile(folder,'simulation_timeshuffle','simulation_data.mat'),fullfile(folder,'simulation_spaceshuffle','simulation_data.mat')};
[sa,X] = network_simulation_beluga.getHeadModel;
location = randi(size(sa.cortex75K.vc,1));
location = randi(size(sa.cortex75K.vc,1),1000,1);
for i = 1:length(folders)
    load(folders{i})
    c1 = corr(dipoles(:,1,1),dipoles(:,1,2));
    c2 = corr(dipoles(:,2,1),dipoles(:,2,2));
    c3 = corr(dipoles(:,3,1),dipoles(:,3,2));
%{
    figureNB
    plot(time,00+squeeze(dipoles(:,1,1)),'color','b'); hold on;

    plot(time,00+squeeze(dipoles(:,1,2)),'color','r')
    plot(time,100+squeeze(dipoles(:,2,1)),'color','b')
    plot(time,100+squeeze(dipoles(:,2,2)),'color','r')
    plot(time,200+squeeze(dipoles(:,3,1)),'color','b')
    plot(time,200+squeeze(dipoles(:,3,2)),'color','r')
    xlim([5e3,10e3])
    text(4900,0+10,'Qx','FontSize',7,'HorizontalAlignment','right')
    text(4900,0-10,sprintf('\\rho=%.2f',c1),'FontSize',7,'HorizontalAlignment','right')
    text(4900,100+10,'Qy','FontSize',7,'HorizontalAlignment','right')
    text(4900,100-10,sprintf('\\rho=%.2f',c2),'FontSize',7,'HorizontalAlignment','right')
    text(4900,200+10,'Qz','FontSize',7,'HorizontalAlignment','right')
    text(4900,200-10,sprintf('\\rho=%.2f',c3),'FontSize',7,'HorizontalAlignment','right')
    axis off;
    drawnow;
%}
    dp = resample(dipoles,1e3,16e3);
    P = [];
    for k = 1:length(location)
        eeg = network_simulation_beluga.getEEG(dp,sa,location(k));
        [psd,f] = pmtm(detrend(eeg,'constant'),2,[],1e3);
        P(:,end+1) = mean(psd,2);
    end

    cDipoles(i,:) = [c1,c2,c3];
    eegPSD(:,i) = mean(P,2);
end

dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';
load(fullfile(dataFolder,'data_time_information.mat'))
t0 = timeInfo.infusion_onset-timeInfo.object_drop;
data = load(fullfile(dataFolder,'data_Cz_multitaper_meanRef.mat'));
pre = [];
for i = 1:length(t0)
    pre(:,i) = nanmedian(data.psd(:,data.time<t0(i),i),2);
end
freq = data.freq;
pre(and(freq>55,freq<65),:) = nan;

load(fullfile(dataFolder,'anatomy_cortical_pairwise_distance_distribution.mat'));
signed_area = A;
total_area = B;
N = 16e9;
dMids = 0.5*(rValues(2:end)+rValues(1:end-1));
nrnCount = mean(diff(signed_area),2)*200000;
nrnCount(end) = N-sum(nrnCount(1:end-1));
corr_kernel = @(d) exp(-d.^2/4);
rho_bar = sum(corr_kernel(dMids).*nrnCount)/sum(nrnCount');
SIG_N = @(rho) N+N*(N-1)*rho*rho_bar;

rho = 0.2;
figureNB(6.1,3.2);
subplot(1,2,1);
    plot(f,smoothdata(eegPSD(:,1),'movmean',10),'color',clrs(1,:));
    hold on;
    plot(f,smoothdata(eegPSD(:,2),'movmean',10),'color',clrs(2,:));
    plot(f,smoothdata(eegPSD(:,3),'movmean',10),'color',clrs(3,:));
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xlabel('Frequency (Hz)');
    xticks([1,10,100])
    xticklabels([1,10,100])
    ylabel('PSD (uV^2/Hz)');
    xax = get(gca,'xaxis');
    xax.TickLabelRotation=0;
    gcaformat;
subplot(1,2,2);
    plotwitherror(freq,pre,'CI','color',[0.5,0.5,0.5]);
    plot(f,smoothdata(eegPSD(:,1),'movmean',10)*SIG_N(0.015),'color',clrs(1,:));
    % plot(f,smoothdata(eegPSD(:,2),'movmean',10)*SIG_N(0),'color',clrs(2,:));
    % plot(f,smoothdata(eegPSD(:,3),'movmean',10)*SIG_N(0),'color',clrs(3,:));
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xlabel('Frequency (Hz)');
    xticks([1,10,100])
    xticklabels([1,10,100])
    ylabel('PSD (uV^2/Hz)');
    xax = get(gca,'xaxis');
    xax.TickLabelRotation=0;
    gcaformat;

figureNB(3.8,2.9);
    dh = 0.1;
    dx = [-dh*2.5,0,dh*2.5];
    for i = 1:3
        for j = 1:3
            mm = cDipoles(i,j);
            fill((i+dx(j))+[-dh,dh,dh,-dh],[0,0,mm,mm],clrs(i,:),'linestyle','none');
            hold on;
        end
    end
    xlim([0.5,3.5])
    ylabel('Dipole correlation')
    gcaformat
    ylim([-0.05,0.6]);
    set(get(gca,'xaxis'),'visible','off');
    line(get(gca,'xlim'),[0,0],'color','k')
    %     B = bar(cDipoles,'FaceColor',[0.5,0.5,0.5],'LineStyle','none');
    % xticklabels({'Constructed input','Shuffled synapses','Poisson'})




for i = 1:length(folders)
    load(folders{i})
    vCorr(i) = corr(V(:,1),V(:,2));
end

figureNB
bar(vCorr)
ylabel('Voltage correlation')
xticklabels({'ST synchrony','T-shuffle','S-shuffle'})
gcaformat






load(folders{1})
figureNB(4.9,4.1);
axes('Position',[0,0,1,1]);
    dp = dipoles;
    dp(time<5e3,:,:) = nan;
    plot(time,00+dp(:,1,2),'color','k'); hold on;
    plot(time,00+dp(:,1,1),'color','r'); hold on;
    plot(time,150+dp(:,2,2),'color','k')
    plot(time,150+dp(:,2,1),'color','r')
    plot(time,300+dp(:,3,2),'color','k')
    plot(time,300+dp(:,3,1),'color','r')

    % xlim([5e3,10e3])
    text(4900,0,'Qx -','FontSize',7,'HorizontalAlignment','right','VerticalAlignment','middle')
    % text(4900,0-10,sprintf('\\rho=%.2f',c1),'FontSize',7,'HorizontalAlignment','right')
    text(4900,150,'Qy -','FontSize',7,'HorizontalAlignment','right','VerticalAlignment','middle')
    % text(4900,100-10,sprintf('\\rho=%.2f',c2),'FontSize',7,'HorizontalAlignment','right')
    text(4900,300,'Qz -','FontSize',7,'HorizontalAlignment','right','VerticalAlignment','middle')
    % text(4900,200-10,sprintf('\\rho=%.2f',c3),'FontSize',7,'HorizontalAlignment','right')
    axis off;
    drawnow;




load('E:\Research_Projects\004_Propofol\data\simulations\raw\synthetic_spikes\synthetic_spikes_analyzed.mat')
R = repmat(1:6,1,11);
for i = 1:length(R)
    C(i,1) = corr(dipoles(:,1,1,i),dipoles(:,1,2,i));
    C(i,2) = corr(dipoles(:,2,1,i),dipoles(:,2,2,i));
    C(i,3) = corr(dipoles(:,3,1,i),dipoles(:,3,2,i));
end


figureNB;
plot(R,nanmean(C,2),'.k','MarkerSize',20)
hold on;
plot(1:6,splitapply(@nanmean,nanmean(C,2)',R),'-k');

cr = reshape(nanmean(C,2),[6,11])';

rRange = 10.^linspace(-4,-1,6);

m0 = nanmean(cr);
s = stderror(cr);
figureNB;
    plot(rRange,m0,'.-k','LineWidth',1,'MarkerSize',10)
    hold on;
    line([rRange;rRange],[m0-1.96*s;m0+1.96*s],'color','k','LineWidth',1)
    set(gca,'xscale','log');
    xlim([1e-4,0.0252])
    ylabel('Dipole correlation')
    xlabel('Synapse correlation (R_{max})')

figureNB(5.25,1)
X = eeg(1.1e5:1.4e5,:);
axes('Position',[0,0,1,1])
    plot(X(:,1),'color','r')
    hold on;
    plot(X(:,2),'color','k')
    line([38e3-250*16,38e3],[-0.5,-0.5]*1e-6,'color','k','LineWidth',1)
    line([38e3-250*16,38e3-250*16],[-0.25,0.75]*1e-6,'color','k','LineWidth',1)
    xlim([1,40e3]);
    axis off;