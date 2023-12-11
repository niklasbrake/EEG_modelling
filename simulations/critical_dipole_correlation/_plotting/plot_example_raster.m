[sa,X] = network_simulation_beluga.getHeadModel;

folder = fullfile(dataFolder,'simulations','example_embedding');
load(fullfile(folder,'model.mat'));
network = network.rebaseNetwork(folder);

[ids,ts,ei] = network.getprenetwork(network.spikingFile);
cons = csvread(fullfile(folder,'postsynaptic_network','connections.csv'));
[~,I] = sort(cons(:,3));
cons = cons(I,:);
idcs = [];
ids2 = []; ts2 = [];
for i = 1:size(cons,1)
    idcs = find(ids==cons(i,3));
    ids2 = [ids2;i+0*idcs];
    ts2 = [ts2;ts(idcs)];
end

mData = network.getmData;

P = arrayfun(@(i,j) mData{i}.pos(j,:),cons(:,1),cons(:,2),'UniformOutput',false);
P = cat(1,P{:});
[thet,phi] = cart2sph(P(:,1),P(:,2),P(:,3));

idx=randi(size(thet,1));
D = network_simulation_beluga.haversine_distance2([phi(idx),thet(idx)],[phi,thet]);
[~,J] = sort(D);
[~,I] = sort(J);
ids3 = I(ids2);


idcs = find(cons(ids3,1)==1);
ids_nrn1 = ids3(idcs);
ts_nrn1 = ts2(idcs);
ids_nrn1 = findgroups(ids_nrn1);

idcs = find(cons(ids3,1)==2);
ids_nrn2 = ids3(idcs);
ts_nrn2 = ts2(idcs);
ids_nrn2 = findgroups(ids_nrn2)+max(ids_nrn1);


load(fullfile(folder,'simulation','simulation_data.mat'));

i0 = 31423;
eeg = network_simulation_beluga.getEEG(dipoles,sa,i0);

figureNB(5.4,3.5);
axes('Position',[0.1,0.77,0.84,0.2]);
    plot(time*1e-3,eeg(:,1)*1e6,'k');
    hold on;
    plot(time*1e-3,eeg(:,2)*1e6,'r');
    xlim([0,4]);
    line([11,12],[2,2],'color','k','LineWidth',1.5)
    axis off;
    yl = get(gca,'ylim');
axes('Position',[0.05,0.77,0.015,0.2]);
    ylim(yl);
    line([0,0],[-1,1],'color','k','linewidth',1.5)
    set(get(gca,'xaxis'),'visible','off');
    axis off;
axes('Position',[0.1000 0.0782 0.8400 0.6418]);
    raster(ids_nrn1,ts_nrn1*1e-3,gcf);
    hold on;
    R = raster(ids_nrn2,ts_nrn2*1e-3,gcf);
    R.Color='r'
    xlim([0,4]);
    xlabel('')
    ylim([1,max(ids_nrn2)]);
    yticks([1,max(ids_nrn1),max(ids_nrn2)]);
    line(get(gca,'xlim'),[1,1]*max(ids_nrn1),'color','k');
    xticks(1:12)
    xticklabels({})
    yticks([])
    ylabel('Synaptic input raster')
    gcaformat;
    set(get(gca,'yaxis'),'visible','on')
    box on;