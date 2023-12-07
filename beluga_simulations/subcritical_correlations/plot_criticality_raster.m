folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\example_embedding';
network = network_simulation_beluga(folder);
network = network.initialize_postsynaptic_network(2,[1,2]);
network.tmax = 50e3;
% load('E:\Research_Projects\004_Propofol\data\simulations\raw\test\connected_components\CC1\model.mat')
N = network.getsynapsecount;

% Get synapse spikes times
cons = csvread(fullfile(network.postNetwork,'connections.csv'));
[ids,ts,ei] = network.getprenetwork(network.spikingFile);
[~,I] = sort(cons(:,3));
cons = cons(I,:);
% nullIDs = find(ids>N);
% ids(nullIDs) = [];
% ts(nullIDs) = [];
for i = 1:30e3
    if(sum(cons(:,3)==i)==0)
        idcs = find(ids==i);
        ids(idcs) = [];
        ts(idcs) = [];
    end
end

% Get synapse locations
mData = network.getmData;
P = arrayfun(@(i,j) mData{i}.pos(j,:),cons(:,1),cons(:,2),'UniformOutput',false);
P = cat(1,P{:});
[thet,phi] = cart2sph(P(:,1),P(:,2),P(:,3));
nrnID = cons(I,1);

idx=randi(size(thet,1));
D = network_simulation_beluga.haversine_distance2([phi(idx),thet(idx)],[phi,thet]);
[~,J] = sort(D);
[~,I] = sort(J);
ids = I(ids);

% raster(ids,ts);

idcs = find(nrnID(ids)==1);
ids1 = ids(idcs);
ts1 = ts(idcs);
ids1 = findgroups(ids1);

idcs = find(nrnID(ids)==2);
ids2 = ids(idcs);
ts2 = ts(idcs);
ids2 = findgroups(ids2)+max(ids1);

% thet1 = thet(unique(ids1));
% phi1 = phi(unique(ids1));
% thet2 = thet(unique(ids2));
% phi2 = phi(unique(ids2));


% Get pdf of input synchronization
h = histcounts(ts,'BinEdges',0.5:network.tmax+0.5);

t = 1:network.tmax;
figureNB(5,4);
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
    raster(ids1,ts1,gcf);
    hold on;
    R = raster(ids2,ts2,gcf);
    R.Color='r'
    xlim([0,5e3]);
    xticklabels([0:5]);
    ylim([1,max(ids2)]);
    yticks([1,max(ids1),max(ids2)]);
    line(get(gca,'xlim'),[1,1]*max(ids1),'color','k');
    xlabel('Time (s)');
    yticks([])
    ylabel('Synaptic input raster')
    gcaformat;
    set(get(gca,'yaxis'),'visible','on')
    box on;