mTypeSegmentationData = fullfile(network_simulation_beluga.resourceFolder,'cortical_column_Hagen','morphology_segmentations.mat');
load(mTypeSegmentationData)
mTypes ={'L23E_oi24rpy1';'L23I_oi38lbc1';'L4E_53rpy1';'L4E_j7_L4stellate';'L4I_oi26rbc1';'L5E_j4a';'L5E_oi15rpy4';'L5I_oi15rbc1';'L6E_51_2a_CNG';'L6E_oi15rpy4';'L6I_oi15rbc1'};
mData = nrnSegs.(mTypes{2});
X1 = [mean(mData.x,2),mean(mData.y,2),mean(mData.z,2)];
mData = nrnSegs.(mTypes{5});
X2 = [mean(mData.x,2),mean(mData.y,2),mean(mData.z,2)];

figureNB;
    [xs,ys,zs] = sphere(25);
    S = surf(xs*0.99,ys*0.99,zs*0.99,0*zs+0.9,'LineStyle',':','FaceAlpha',0.7,'EdgeColor',[0.6,0.6,0.6]);
    hold on;
    plotProjectDendrites(2,'k');
    plotProjectDendrites(5,'r');
    view([0,0]);
    set(gca,'DataAspectRatio',[1,1,1]);
    gcaformat
    axis off;
    set(gca,'CLim',[0,1]);
    colormap('gray');
return;
fig1 = render_neuron_morphology(2)
xl1 = get(gca,'xlim');
yl1 = get(gca,'ylim');
fig2 = render_neuron_morphology(5)
xl2 = get(gca,'xlim');
yl2 = get(gca,'ylim');

xl = [min(xl1(1),xl2(1)),max(xl1(2),xl2(2))];
yl = [min(yl1(1),yl2(1)),max(yl1(2),yl2(2))];
figure(fig1);
xlim(xl);
ylim(yl);

figure(fig2);
xlim(xl);
ylim(yl);


colormap([1,0,0])

function plotProjectDendrites(iM,clr)
    mTypes ={'L23E_oi24rpy1';'L23I_oi38lbc1';'L4E_53rpy1';'L4E_j7_L4stellate';'L4I_oi26rbc1';'L5E_j4a';'L5E_oi15rpy4';'L5I_oi15rbc1';'L6E_51_2a_CNG';'L6E_oi15rpy4';'L6I_oi15rbc1'};
    mType = mTypes{iM};

    folder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data\cortical_column_Hagen\matlab_morph_data';
    load(fullfile(folder,[mType '_connections.mat']));
    connections = connections+1;

    X = csvread(fullfile(folder,[mType '_segments.csv']));
    for i=  1:size(X,1)
        segs{i} = X(i,X(i,:)~=0);
    end

    xSoma = data(data(:,2)==1,3:5);
    data(:,3:5) = data(:,3:5) - mean(xSoma);
    data(:,6) = max(data(:,6),0.25);

    for i = 1:length(segs)
        idcs = segs{i};
        xSegment = data(idcs,3:6);
        somaIdcs = find(data(idcs,2)==1);
        xSegment(somaIdcs,:) = [];
        lw = max(xSegment(:,4));
        xSegment = xSegment./vecnorm(xSegment,2,2);
        plot3(xSegment(:,1),xSegment(:,2),xSegment(:,3),'.','color',clr,'MarkerSize',7)
    end
end
