function optimize_synapse_geometry(dist_mat)

M = 500;
N = M/2;

subIdcs = randperm(M,N);
dist_mat = dist_mat(subIdcs,:);
dist_mat = dist_mat(:,subIdcs);

flat_dist_mat = squareform(dist_mat);
res_linkage = linkage(flat_dist_mat);
res_order = seriation(res_linkage,N,2*N-1);

seriated_dist = dist_mat;
seriated_dist = seriated_dist(res_order,:);
seriated_dist = seriated_dist(:,res_order);

figureNB(20,6.6);
ax = subplot(1,2,1); hold on
    ax.ColorOrder = clrsPT.lines(10);
    plot(Y(:,1:10)+[1:10]*5)
    yticks([1:10]*5)
    yticklabels([1:10])
    axis ij
    ylim([0,55])
    ylabel('Neuron count')
    xlabel('Time');
    gcaformat_dark;
subplot(1,2,2);
    imagesc((1-seriated_dist).*(1-eye(N)));
    set(gca,'CLim',[0,0.25])
    axis square;
    colormap(clrsPT.sequential(1e3))
    C = colorbar;
    C.Label.String = 'Pairwise correlation'
    C.Color = 'w'
    xlabel('Neuron count')
    ylabel('Neuron count')
    gcaformat_dark;

file = 'E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\network_embedding\file.csv';
csvwrite(file,dist_mat);
pyFun = 'E:\Research_Projects\004_Propofol\code\simulations\functions\embed_data.py';
[err,prints] = system(['python "' pyFun '" ' file]);

data = csvread(file);
x = cos(data(:,1)).*sin(data(:,2));
y = sin(data(:,1)).*sin(data(:,2));
z = cos(data(:,2));
embedding = [x(:),y(:),z(:)]';

leftIdcs = setdiff(1:M,subIdcs);
X = zeros(3,M);
X(:,subIdcs) = embedding;
for i = 1:length(leftIdcs)
    w = CC(leftIdcs(i),subIdcs);
    w(w<quantile(w,0.9)) = 0;
    X(:,leftIdcs(i)) = wgeodmeanSm(embedding,w);
end
theta2 = atan2(X(1,:)',X(2,:)');
phi2 = -acos(X(3,:)');

mTypeSegmentationData = 'E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\cortical_column_Hagen\segment_areas.mat';
load(mTypeSegmentationData)
fn = fieldnames(nrnSegSA);
X = []; Y = [];  Z = []; ID = [];
for i = 1%:length(fn)
    X = [X;nrnSegSA.(fn{i}).vec2soma(:,1)];
    Y = [Y;nrnSegSA.(fn{i}).vec2soma(:,2)];
    Z = [Z;nrnSegSA.(fn{i}).vec2soma(:,3)];
    ID = [ID;i*ones(size(nrnSegSA.(fn{i}).vec2soma(:,3)))];
end
dends = [X,Y,Z];
theta = atan2(X,Y);
phi = -acos(Z);

segs = randperm(size(dends,1),M);
preSyns = ones(M,1);
t2 = data(:,1);
p2 = data(:,2);
syn = [];
for i = 1:length(segs)
    t1 = theta(segs(i));
    p1 = phi(segs(i));
    D = 1-haversine_distance([t1,p1],[theta2,phi2]);
    [~,syn(i)] = max(D.*preSyns);
    preSyns(syn(i)) = 0;
end

clrs = clrsPT.lines(m);
fig = figureNB;
    [xs,ys,zs] = sphere(25);
    sph = mesh(xs,ys,zs); hold on;
    sph.LineStyle = ':'
    sph.EdgeColor = [1,1,1]*0.6;
    sph.FaceColor = [1,1,1]*0.1;
    sph.FaceAlpha = 0.5;
for i = 1:m
    j = segs(find(sum(syn(:)==idcs(i,:),2)));
    plot3(cos(theta(j)).*sin(phi(j)),sin(theta(j)).*sin(phi(j)),cos(phi(j)),'.','color',clrs(i,:),'MarkerSize',15);
    hold on
end
set(gca,'DataAspectRatio',[1,1,1]);
axis off;
gcaformat_dark;

clrs = clrsPT.lines(m);
fig = figureNB;
    [xs,ys,zs] = sphere(25);
    sph = mesh(xs,ys,zs); hold on;
    sph.LineStyle = ':'
    sph.EdgeColor = [1,1,1]*0.6;
    sph.FaceColor = [1,1,1]*0.1;
    sph.FaceAlpha = 0.5;
for i = 1:m
    j = idcs(i,:);
    plot3(cos(theta2(j)).*sin(phi2(j)),sin(theta2(j)).*sin(phi2(j)),cos(phi2(j)),'.','color',clrs(i,:),'MarkerSize',15);
    hold on
end
set(gca,'DataAspectRatio',[1,1,1]);
axis off;
gcaformat_dark;

figureNB;
plot3(nrn.pos(:,1),nrn.pos(:,2),nrn.pos(:,3),'.w');
hold on;
for i = 1:m
    j = segs(find(sum(syn(:)==idcs(i,:),2)));
    plot3(nrn.pos(j,1),nrn.pos(j,2),nrn.pos(j,3),'.','color',clrs(i,:),'MarkerSize',15);
    hold on 
end
set(gca,'DataAspectRatio',[1,1,1]);
axis off;
gcaformat_dark;

function d = hav(x)
    d = (1-cos(x))/2;

end
function d = haversine_distance(p1,x)
    d = hav(p1(1)-x(:,1))+(1-hav(x(:,1)-p1(1))-hav(p1(1)+x(:,1))).*hav(p1(2)-x(:,2));
end

end