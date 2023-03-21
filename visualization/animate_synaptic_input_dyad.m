masterPath = 'E:\Research_Projects\004_Propofol\data\simulations\raw\test';
% masterPath = 'E:\Research_Projects\004_Propofol\data\simulations\raw\multisynapse\m=0.98';


fid = fopen(fullfile(masterPath,'postsynaptic_network','mTypes.txt'));
mType2 = textscan(fid,'%s');
mType2 = mType2{1};
fclose(fid);
for i = 1:length(mType2)
    [~,mType,~] = fileparts(mType2{i});
    load(fullfile(network_simulation_beluga.resourceFolder,'cortical_column_Hagen','segment_areas.mat'));
    mData{i} = nrnSegs.(mType);
    mData{i}.x = mData{i}.x+(i-1)*500;
    mData{i}.pos = [mean(mData{i}.x,2),mean(mData{i}.y,2),mean(mData{i}.z,2)];
end

cons = csvread(fullfile(masterPath,'postsynaptic_network','connections.csv'));


[ids,ts,ei] = network_simulation_beluga.getprenetwork(fullfile(masterPath,'presynaptic_network','spikeTimes.csv'));

N = length(ei);

loc = nan(N,3);
for i = 1:N
    seg = find(cons(:,3)==i);
    if(~isempty(seg))
        seg = cons(seg,2);
        nrnID = cons(seg,1);
        loc(i,:) = mData{nrnID}.pos(seg,:);
    end
end

loc2 = csvread(fullfile(masterPath,'presynaptic_network','locations.csv'));
umap = network_simulation_beluga.loadSynapseLocations(fullfile(masterPath,'presynaptic_network','UMAP_embedding.csv'),N);
x = umap(:,1);
y = umap(:,2);
z = umap(:,3);

parents = csvread(fullfile(masterPath,'presynaptic_network','multisynapse_IDs.csv'));
idcs = find(parents~=-1);
% loc2(idcs,:) = nan;

% figureNB;
%     S2 = scatter(loc2(:,1),loc2(:,2),ones(N,1),nan(N,1),'filled');
%     set(gca,'CLim',[0,1]);
%     set(gca,'DataAspectRatio',[1,1,1]);
%     gcaformat_dark;
%     axis off;
%     colormap([1,1,0;0.5,0.5,0.5;0,0,1]);
%     for t = 0:4:max(ts)
%         idcs = ids(find(and(ts>t,ts<t+4)));
%         S2.CData(idcs) = ei(idcs);
%         S2.SizeData(idcs) = 8;
%         drawnow;
%         pause(max(0.05-toc,0));
%         S2.CData(idcs) = nan;
%         S2.SizeData(idcs) = 1;
%         tic
%     end


figureNB(24,9.7);
ax1 = subplot(1,4,1);
    S2 = scatter(loc2(:,1),loc2(:,2),ones(N,1),0.5*ones(N,1),'filled');
    set(gca,'CLim',[0,1]);
    set(gca,'DataAspectRatio',[1,1,1]);
    gcaformat_dark;
    axis off;
    colormap(ax1,[1,1,0;0.5,0.5,0.5;0,0,1]);
ax2 = subplot(1,4,2:4);
    S = scatter3(loc(:,1),loc(:,2),loc(:,3),3*ones(N,1),nan(N,1),'filled');
    set(gca,'CLim',[0,1]);
    set(gca,'DataAspectRatio',[1,1,1]);
    for i = 1:length(mData)
        line(mData{i}.x',mData{i}.y',mData{i}.z','color','w');
    end
    gcaformat_dark;
    axis off;
    view([0,0])
    colormap(ax2,[1,1,0;,0,0,1]);

tic;
for t = 0:4:max(ts)
    idcs = ids(find(and(ts>t,ts<t+4)));
    S.CData(idcs) = ei(idcs);
    S2.CData(idcs) = ei(idcs);
    S.SizeData(idcs) = 8;
    S2.SizeData(idcs) = 8;
    drawnow;
    pause(max(0.05-toc,0));
    S.CData(idcs) = nan;
    S2.CData(idcs) = 0.5;
    S.SizeData(idcs) = 1;
    S2.SizeData(idcs) = 1;
    tic
end

