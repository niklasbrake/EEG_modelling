folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\example_embedding'
[sa,X] = network_simulation_beluga.getHeadModel;

load(fullfile(folder,'model.mat'));
[ids,ts,ei] = network.getprenetwork(network.spikingFile);

X = csvread(fullfile(folder,'presynaptic_network\locations.csv'));


cons = csvread(fullfile(network.postNetwork,'connections.csv'));
nrnID = cons(:,1);

nrnID = nan(30e3,1);
nrnID(cons(:,3)) = cons(:,1);


[~,I] = sort(X(:,1));
[~,J] = sort(I);


idcs = find(nrnID(ids)==1);
ids1 = J(ids(idcs));
ts1 = ts(idcs);
ids1 = findgroups(ids1);

idcs = find(nrnID(ids)==2);
ids2 = J(ids(idcs));
ts2 = ts(idcs);
ids2 = findgroups(ids2)+max(ids1);


load(fullfile(folder,'simulation\simulation_data.mat'))

eeg = network_simulation_beluga.getEEG(dipoles,sa,7e3);
corr(eeg)

figureNB(5.4,3.5);
axes('Position',[0.16, 0.77, 0.77, 0.19])
    plot(time,eeg(:,1),'color','k');
    hold on;
    plot(time,eeg(:,2),'color','r');
    xlim([4e3,6e3])
    % ylim([-0.5,0.5]*1e-6);
    ylim([-0.8,0.8]*1e-6)
    axis off;
    yh = 0.19/1.6;
    annotation('line','Position',[0.13,0.77+0.19/2-yh/2,0,yh],'LineWidth',0.75)
    line([5.75,6]*1e3,[0.8,0.8]*1e-6,'LineWidth',0.75,'color','k');

axes('Position',[0.16, 0.11, 0.77, 0.6])
    R = raster(ids1,ts1,gcf);
    hold on
    R = raster(ids2,ts2,gcf);
    R.Color = 'r';
    ylim([0,max(ids2)])
    xlim([4e3,6e3])
    line([1e3,6e3],[max(ids1),max(ids1)],'color','k','LineWidth',0.75)
    xticks([2e3:250:6e3]);
    yticks([]);
    set(get(gca,'yaxis'),'visible','on')
    set(get(gca,'xaxis'),'visible','on')
    xlabel('')
    xticklabels({})
    gcaformat;
    box on

