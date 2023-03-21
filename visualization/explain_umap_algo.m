folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\test\correlations';

load('E:\Research_Projects\004_Propofol\data\resources\cortical_column_Hagen\segment_areas.mat');
nrnID = 1;
fn = fieldnames(nrnSegs);
mData = nrnSegs.(fn{nrnID});

N = 500;
% [dist_list,~,k] = getToyModel(N);
network = network_simulation(folder);
network = network.initialize_postsynaptic_network(1,nrnID);
network = network.setsynapsecount(N);
% csvwrite(fullfile(network.preNetwork,'correlations.csv'),dist_list);
network.form_connections(1);

synCons = csvread(fullfile(network.postNetwork,'connections',[network.neurons(1).name '.csv']));

N0 = 200;
[~,seriated_dist] = getToyModel(N0);


blue = clrsPT.qualitative_CM.blue;
red = clrsPT.qualitative_CM.orange;
clrs0 = [blue;red];
figureNB(11.5,3);
ax(1) = axes('Position',[0.07,0.27,0.14,0.55]);
    imagesc((1-seriated_dist).*(1-eye(N0))); hold on;
        axis square;
        axis xy;
        xlim([0,N0])
        ylim([0,N0])
    F = fill([0,N0,N0,0],[0,0,N0,N0],'w','LineWidth',0.75);
        F.EdgeColor = 'k';
        F.FaceAlpha = 0;
    F = fill([0,N0/2,N0/2,0],[0,0,N0/2,N0/2],'w','LineWidth',1.5);
        F.EdgeColor = clrs0(1,:);
        F.FaceAlpha = 0;
    F = fill([N0/2+1,N0,N0,N0/2+1],[N0/2+1,N0/2+1,N0,N0],'w','LineWidth',1.5);
        F.EdgeColor = clrs0(2,:);
        F.FaceAlpha = 0;
    C = colorbar('Position',[0.23,0.27,0.02,0.55]);;
        C.Color = 'k';
        C.Ticks = [];
        CM = gray(1e3);
        colormap(flip(CM(200:end,:)));
        set(gca,'CLim',[0,0.25])
    xlabel('Neuron count')
    ylabel('Neuron count')
    axis off;
    gcaformat;
    text(N0/2+1,-5,'Neuron count','FontSize',6,'HorizontalAlignment','center','VerticalAlignment','top','color','k')
    text(-10,N0/2+1,'Neuron count','FontSize',6,'HorizontalAlignment','center','VerticalAlignment','bottom','rotation',90,'color','k')
    gcaformat;
ax(2)=axes('Position',[0.26,0.03,0.4,1]);
    idcs = randperm(500,200);
    [xs,ys,zs] = sphere(25);
    ms = mesh(xs,ys,zs); hold on;
    ms.LineStyle = ':';
    ms.EdgeColor = 1-[1,1,1]*0.6;
    ms.FaceColor = 1-[1,1,1]*0.1;
    ms.FaceAlpha = 0.5;
    scatter3(umap(idcs,1),umap(idcs,2),umap(idcs,3),5,clrs(idcs,:),'filled');
    set(gca,'DataAspectRatio',[1,1,1]);
    axis off;
ax(3)=axes('Position',[0.62,0.1,0.2,0.9]);
    scatter3(synPos(:,1),synPos(:,2),synPos(:,3),2,clrs(synID,:),'filled'); hold on;
    set(gca,'DataAspectRatio',[1,1,1]);
    line(dendPos.x',dendPos.y',dendPos.z','color','k','LineWidth',0.25);
    xlim([-200,200]);
    ylim([-200,200]);
    zlim([-200,350]);
    view([-90,10]);
    axis off;
ax(4)=axes('Position',[0.8,0.1,0.2,0.9]);
    scatter3(synPos2(:,1),synPos2(:,2),synPos2(:,3),2,clrs(synID2,:),'filled'); hold on;
    set(gca,'DataAspectRatio',[1,1,1]);
    line(dendPos2.x',dendPos2.y',dendPos2.z','color','k','LineWidth',0.25);
    view([-120,10]);
    xlim([-200,200]);
    ylim([-200,200]);
    zlim([-200,350]);
    axis off;


label = @(x,y,str) annotation('textbox', [x,y,0.3,0.05],'String',str, 'LineStyle', 'none', 'FontSize',7,'Margin',0,'HorizontalAlignment','center');

txt = label(0,0.05,'Pairwise correlation');
txt = label(0.31,0.05,'UMAP embedding');
txt = label(0.68,0.05,'Projection onto dendrites');


function [synPos,dendPos,synID,bestD] = getSynapseEmbedding(data,nrnID,k,umap)
    N = size(data,1);
    M = floor(N);
    load('E:\Research_Projects\004_Propofol\data\resources\cortical_column_Hagen\segment_areas.mat');
    mData = nrnSegs.(nrnID);
    pos = [mean(mData.x,2),mean(mData.y,2),mean(mData.z,2)];

    sa = mData.area;
    X = pos(:,1);
    Y = pos(:,2);
    Z = pos(:,3);
    [theta,phi] = cart2sph(X,Y,Z);
    dendriteEmbedding = [phi(:),theta(:)];
    synapseEmbedding = [data(:,3)+pi/2,data(:,2)];
    saCDF = cumsum(sa)/sum(sa);
    iSegs = interp1(saCDF,1:length(saCDF),rand(M,1),'next','extrap');
    % iSegs = randi(length(saCDF),M,1);
    syn = zeros(size(iSegs));
    remainingSyns = true(N,1);
    for i = 1:M
        D = 1-network_simulation.haversine_distance2(dendriteEmbedding(iSegs(i),:),synapseEmbedding);
        [bestD(i),syn(i)] = max(D.*remainingSyns);
        remainingSyns(syn(i)) = false;
    end

    r = rand(size(iSegs));
    synPos(:,1) = mData.x(iSegs,1).*r+mData.x(iSegs,2).*(1-r);
    synPos(:,2) = mData.y(iSegs,1).*r+mData.y(iSegs,2).*(1-r);
    synPos(:,3) = mData.z(iSegs,1).*r+mData.z(iSegs,2).*(1-r);
    synPos = synPos;% + randn(size(synPos));
    dendPos = mData;
    synID = syn;
end

function [dist_list,seriated_dist,k] = getToyModel(N)
    M = 2;
    T = 500;
    drivers = randn(T,M);
    % k = randi(M,N,1);
    k = zeros(N,1);
    idcs = randperm(N,N/2);
    k(idcs) = 1;
    X = zeros(T,N);
    rho = 0.4;
    for i = 1:N
        if(k(i)==1)
            % w(1) = betarnd(5,2);
            w(1) = 1;
            w(2) = 1-w(1);
            w = w*rho;
        else
            % w(1) = betarnd(2,5);
            w(1) = 0;
            w(2) = 1-w(1);
            w = w*rho;
        end
        X(:,i) = sum(w.*drivers,2)+sqrt(1-rho^2)*randn(T,1);
        W(i) = w(1)/rho;
    end
    dist_mat = 1-abs(corr(X));

    flat_dist_mat = squareform(dist_mat);
    res_linkage = linkage(flat_dist_mat);
    res_order = seriation(res_linkage,N,2*N-1);
    seriated_dist = dist_mat;
    seriated_dist = seriated_dist(res_order,:);
    seriated_dist = seriated_dist(:,res_order);


    dist_list = zeros(N*N,3);
    l=1;
    for i = 1:N
        for j = 1:N
            dist_list(l,:) = [i-1,j-1,dist_mat(i,j)];
            l = l+1;
        end
    end
end