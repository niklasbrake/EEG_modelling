folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\critical_embedding';
network = runSimulation(folder);

[synapseEmbedding,idPreSyn] = network_simulation_beluga.loadSynapseLocations(fullfile(network.preNetwork,'UMAP_embedding.csv'),network.getsynapsecount);

locations2D = csvread(fullfile(network.preNetwork,'locations.csv'));
correlations = csvread(fullfile(network.preNetwork,'correlations.csv'));
[ids,ts,ei] = network.getprenetwork(network.spikingFile);

dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';
load(fullfile(dataFolder,'cortical_column_Hagen','morphology_segmentations.mat'));
fn = fieldnames(nrnSegs);
mData = nrnSegs.(fn{1});
mData2 = nrnSegs.(fn{2});

synCons = csvread(fullfile(network.postNetwork,'connections.csv'));
idcs1 = find(synCons(:,1)==1);
r = rand(size(idcs1,1),1);
synPos = [];
synPos(:,1) = mData.x(synCons(idcs1,2),1).*r+mData.x(synCons(idcs1,2),2).*(1-r);
synPos(:,2) = mData.y(synCons(idcs1,2),1).*r+mData.y(synCons(idcs1,2),2).*(1-r);
synPos(:,3) = mData.z(synCons(idcs1,2),1).*r+mData.z(synCons(idcs1,2),2).*(1-r);
synID = synCons(idcs1,3);

idcs2 = find(synCons(:,1)==2);
r = rand(size(idcs2,1),1);
synPos2 = [];
synPos2(:,1) = mData2.x(synCons(idcs2,2),1).*r+mData2.x(synCons(idcs2,2),2).*(1-r);
synPos2(:,2) = mData2.y(synCons(idcs2,2),1).*r+mData2.y(synCons(idcs2,2),2).*(1-r);
synPos2(:,3) = mData2.z(synCons(idcs2,2),1).*r+mData2.z(synCons(idcs2,2),2).*(1-r);
synID2 = synCons(idcs2,3);

[xs,ys,zs] = sphere(25);
[~,I] = sort(locations2D(idPreSyn,2));
[~,I1] = sort(locations2D(synID,2));
[~,I2] = sort(locations2D(synID2,2));

iMax = randi(size(correlations,1));
% iMax = 307223;
J1 = correlations(iMax,1)+1;
J2 = correlations(iMax,2)+1;
idcs1 = find(ids==J1);
idcs2 = find(ids==J2);
idcs = randperm(network.getsynapsecount,network.getsynapsecount-1e4);
idcs3 = [];
idcs4 = [];
for i = 1:length(idcs)
    tempI = find(I1==idcs(i));
    if(~isempty(tempI))
        idcs3 = [idcs3;tempI];
    end
    tempI = find(I2==idcs(i));
    if(~isempty(tempI))
        idcs4 = [idcs4;tempI];
    end
end

CM = hsv(network.getsynapsecount);
CM = 0.3*[255,254,222]/255+0.5*CM;

CM = hot(network.getsynapsecount*1.4);
CM = CM(0.05*network.getsynapsecount+1:1.05*network.getsynapsecount,:);
% CM = clrsPT.sequential(network.getsynapsecount*1.4);
% CM = flipud(CM(end-network.getsynapsecount+1:end,:));

figureNB(11.5,3);
ax(1) = axes('Position',[0.095,0.22,0.2,0.55]);
    line((ts(idcs1).*[1,1])',(ts(idcs1)*0+[0,1])','color',CM(J1,:));
    line((ts(idcs2).*[1,1])',(ts(idcs2)*0+[1,2])','color',CM(J2,:));
    text(-2.2e4,0.5,sprintf('ID %d',J1),'FontSize',6,'Color',CM(J1,:))
    text(-2.2e4,1.5,sprintf('ID %d',J2),'FontSize',6,'Color',CM(J2,:))
    % title(sprintf('STTC = %.3f',C(iMax,3)),'FontSize',7,'FontWeight','normal')
    axis off;
ax(2)=axes('Position',[0.27,0.03,0.4,1]);
    surf(xs,ys,zs,zs*0,'LineStyle',':','FaceAlpha',0.1,'EdgeColor','none')
    hold on;
    S = scatter3(synapseEmbedding(I,1),synapseEmbedding(I,2),synapseEmbedding(I,3),1,CM,'filled');
        S.XData(idcs)  = [];
        S.YData(idcs)  = [];
        S.ZData(idcs)  = [];
        S.CData(idcs,:)  = [];
    set(gca,'DataAspectRatio',[1,1,1]);
    set(gca,'CLim',[0,1]);
    colormap('gray');
    axis off;
    axis vis3d
    view([-120,10]);
ax(3)=axes('Position',[0.62,0.1,0.2,0.9]);
    S = scatter3(synPos(I1,1),synPos(I1,2),synPos(I1,3),1,CM(sort(synID),:),'filled'); hold on;
        S.XData(idcs3)  = [];
        S.YData(idcs3)  = [];
        S.ZData(idcs3)  = [];
        S.CData(idcs3,:)  = [];
    set(gca,'DataAspectRatio',[1,1,1]);
    line(mData.x',mData.y',mData.z','color','k','LineWidth',0.25);
    xlim([-200,200]);
    ylim([-200,200]);
    zlim([-200,350]);
    view([-90,10]);
    axis off;
ax(4)=axes('Position',[0.8,0.1,0.2,0.9]);
    S = scatter3(synPos2(I2,1),synPos2(I2,2),synPos2(I2,3),1,CM(sort(synID2),:),'filled'); hold on;
        S.XData(idcs4)  = [];
        S.YData(idcs4)  = [];
        S.ZData(idcs4)  = [];
        S.CData(idcs4,:)  = [];
    set(gca,'DataAspectRatio',[1,1,1]);
    line(mData2.x',mData2.y',mData2.z','color','k','LineWidth',0.25);
    view([-120,10]);
    xlim([-200,200]);
    ylim([-200,200]);
    zlim([-200,350]);
    axis off;


function network = runSimulation(folder)
    % Initialize network
    network = network_simulation_beluga(folder);

    % Initialize post network
    nrnID = [1,2];
    network = network.initialize_postsynaptic_network(2,nrnID);

    % Presyanptic network parameters
    network.tmax = 50e3; % 2 seconds
    network.branchNo = 0.98;

    if(exist(fullfile(network.preNetwork,'UMAP_embedding.csv')))
        disp('Simulation already run. Returning results.');
        return;
    end
    network.simulatespikes_critplane(network.getsynapsecount,network.tmax);

    % Compute presynaptic correlations
    network.compute_presynaptic_correlations;

    % Perform UMAP
    network.embed_presyanptic_neurons;

    % Optimize synapse placement
    network = network.form_connections(1);
end