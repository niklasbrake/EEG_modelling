clrs = [clrsPT.qualitative_CM.blue;clrsPT.qualitative_CM.red];
clrs = [[112,58,160]/255;[192,0,0]/255];

repCount = 10;
tmax0 = 2e3;
tmax = 2e3;
nrnCount = 2;

N = 1e4;
network = network_simulation_beluga('C:/Users/brake/Documents/temp/network_osc');
network = network.initialize_postsynaptic_network(2);
network.tmax = tmax;
network.branchNo = 0.5;

[ids,ts] = simulatespikes_osc(N,tmax0,network);

idcs = find(ids<=N/2);
ids1 = ids(idcs);
ids0 = randperm(N/2);
ids1 = ids0(ids1);
ts1 = ts(idcs);

idcs = find(ids>N/2);
ids2 = ids(idcs)-N/2;
ids0 = randperm(N/2);
ids2 = ids0(ids2);
ts2 = ts(idcs);


fig = figureNB(4.25,4.25);

[h,bins] = histcounts(ts1,'BinWidth',10); 
h = h/N/10e-3;
t0 = bins(1:end-1)+2;
axes('Position',[0.15,0.5,0.8,0.4]);
    r = raster(ids1,ts1,fig);
    r.Color = clrs(2,:);
    ylim(range(ids1)*[-0.1,1.2]);
    axis off
    line([9e3,10e3],-0.1*max(ids1)*[1,1],'color','k','LineWidth',1)
    text(9500,-0.12*max(ids1),'1 s','FontSize',6,'VerticalAlignment','top','HorizontalAlignment','center');
    xlim([0,750]);
    text(4.9e3,max(ids1)/2,{'Presynaptic','neurons'},'Rotation',90,'FontSize',6,'HorizontalAlignment','center','VerticalAlignment','bottom');
axes('Position',[0.15,0.85,0.8,0.1]);
    plot(t0,h,'color',clrs(2,:),'LineWidth',1);
    xlim([0,750]);
    ylim([0,2])
    yticklabels([]);
    yticks([0,2]);
    gcaformat;
    set(get(gca,'xaxis'),'visible','off');
    text(4.9e3,2.5,{'Hz',''},'Rotation',90,'FontSize',6,'HorizontalAlignment','center','VerticalAlignment','bottom');
    text(4.9e3,0,'0 ','FontSize',6,'HorizontalAlignment','right','VerticalAlignment','middle');
    text(4.9e3,5,'5 ','FontSize',6,'HorizontalAlignment','right','VerticalAlignment','middle');
    axis off;

[h,bins] = histcounts(ts2,'BinWidth',10); 
h = h/N/10e-3;
t0 = bins(1:end-1)+2;
axes('Position',[0.15,0.05,0.8,0.4]);
    r = raster(ids2,ts2,fig);
    r.Color = clrs(1,:);
    ylim(range(ids2)*[-0.1,1.2]);
    axis off
    line([9e3,10e3],-0.1*max(ids2)*[1,1],'color','k','LineWidth',1)
    text(9500,-0.12*max(ids2),'1 s','FontSize',6,'VerticalAlignment','top','HorizontalAlignment','center');
    xlim([0,750]);
    text(4.9e3,max(ids2)/2,{'Presynaptic','neurons'},'Rotation',90,'FontSize',6,'HorizontalAlignment','center','VerticalAlignment','bottom');
axes('Position',[0.15,0.4,0.8,0.1]);
    plot(t0,h,'color',clrs(1,:),'LineWidth',1);
    xlim([0,750]);
    ylim([0,2])
    yticklabels([]);
    yticks([0,2]);
    gcaformat;
    set(get(gca,'xaxis'),'visible','off');
    text(4.9e3,2.5,{'Hz',''},'Rotation',90,'FontSize',6,'HorizontalAlignment','center','VerticalAlignment','bottom');
    text(4.9e3,0,'0 ','FontSize',6,'HorizontalAlignment','right','VerticalAlignment','middle');
    text(4.9e3,5,'5 ','FontSize',6,'HorizontalAlignment','right','VerticalAlignment','middle');
    axis off;




load(fullfile(network_simulation_beluga.resourceFolder,'cortical_column_Hagen','segment_areas.mat'))
mData = nrnSegs.L23E_oi24rpy1;
synPos = csvread('C:\Users\brake\Documents\temp\network_osc\postsynaptic_network\UMAP_embedding.csv');
synData = csvread('C:\Users\brake\Documents\temp\network_osc\postsynaptic_network\connections.csv');
synData = synData(synData(:,1)==1,:);

segID = synData(:,2);
clr = synPost(synData(:,3),3)>=0;
idcs = randperm(10717);
segID = segID(idcs);
clr = clr(idcs,:);
pos = [mean(mData.x,2),mean(mData.y,2),mean(mData.z,2)];

figureNB(3,4);
axes('Position',[0,0,1,1]);
    line(mData.x',mData.z','color','k','LineWidth',0.5);
    hold on;
    scatter(pos(segID,1),pos(segID,3),0.5,clr,'filled');
    colormap(clrs);
    xticks([])
    yticks([])
    set(gca,'DataAspectRatio',[1,1,1]);
    axis off;



load('C:\Users\brake\Documents\temp\network_osc\postsynaptic_network\simulation_data.mat');
t = time(2:end);
y = dipoles(2:end,:,2);
y(t<1e3,:) = nan;

offset = linspace(0,50,3);
figureNB(4,4.5);
    plot(time(2:end),y+offset,'color','k')
    xlim([800,2e3])
    
    line([900,950],[1,1]*offset(1),'color','k','LineWidth',1);
    text(870,offset(1),'Px','FontSize',6,'color','k','HorizontalAlignment','right');
    
    line([900,950],[1,1]*offset(2),'color','k','LineWidth',1);
    text(870,offset(2),'Py','FontSize',6,'color','k','HorizontalAlignment','right');

    line([900,950],[1,1]*offset(3),'color','k','LineWidth',1);
    text(870,offset(3),'Pz','FontSize',6,'color','k','HorizontalAlignment','right');

    line([1500,1700],[-20,-20],'LineWidth',1.5,'color','k');
    text(1600,-22,'200 ms','FontSize',6,'HorizontalAlignment','center','VerticalAlignment','top')

    line([1400,1400],[-20,-10],'LineWidth',1.5,'color','k');
    text(1350,-15,['10 nA' char(956) 'm'],'FontSize',6,'HorizontalAlignment','right','VerticalAlignment','middle')






tauResults = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\analzye_simulations_tau_osc.mat');
tauResults.P1 = tauResults.P1(:,:,[1,end]);
mResults = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\analzye_simulations_crit_osc.mat');
mResults.P1 = mResults.P1(:,:,[1,5]);
taumResults = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\analyze_simulations_mixed_osc_incomplete.mat');
taumResults.P1 = taumResults.P;
f = mResults.f;
% clrs = [clrsPT.qualitative_CM.blue;clrsPT.qualitative_CM.red];
clrs = [0,0,0;clrsPT.qualitative_CM.red];
% clrs = clrs([1,end],:);

m = [0,0.98];
figureNB(11.4,4);
axes('Position',[0.11,0.28,0.19,0.65]);
    h=[];
    for i = 1:length(m)
        y = mResults.P1(:,:,i);
        % y = y./y(f==100,:);
        h(i) = plotwitherror(f,y,'CI','color',clrs(i,:),'LineWidth',1);
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylim(10.^[-17.2,-13.5])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    L = legend(h,{'m = 0','m = 0.98'}); L.Box = 'off'; L.ItemTokenSize = [10,5]; 
    gcaformat;
    % T = title({'Network criticality','(\tau = 10 ms)'},'fontsize',6,'FontWeight','normal');
    % T.Position(2) = 3e-14;
    L.Position(1:2) = [0.11,0.28];
tau = [10,35];
axes('Position',[0.44,0.28,0.19,0.65]);
    h=[];
    for i = 1:length(tau)
        h(i) = plotwitherror(f,tauResults.P1(:,:,i),'CI','color',clrs(i,:),'LineWidth',1);
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylim(10.^[-17.2,-13.5])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    L = legend(h,{'\tau = 10 ms','\tau = 30 ms'}); L.Box = 'off'; L.ItemTokenSize = [10,5]; 
    gcaformat;
    % T = title({'Synaptic decay','(m=0)'},'fontsize',6,'FontWeight','normal');
    % T.Position(2) = 3e-14;
    L.Position(1:2) = [0.44,0.28];
axes('Position',[0.77,0.28,0.19,0.65]);
    h=[];
    for i = 1:length(tau)
        h(i) = plotwitherror(f,taumResults.P1(:,:,i),'CI','color',clrs(i,:),'LineWidth',1);
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylim(10.^[-17.2,-13.5])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    gcaformat;
    L = legend(h,{'\tau = 10 ms','\tau = 30 ms'}); L.Box = 'off'; L.ItemTokenSize = [10,5]; 
    gcaformat;
    % T = title({'Synaptic decay','(mixed network)'},'fontsize',6,'FontWeight','normal');
    % T.Position(2) = 3e-14;
    L.Position(1:2) = [0.77,0.28];