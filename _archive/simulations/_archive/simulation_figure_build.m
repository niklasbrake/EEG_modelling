

load('E:\Research_Projects\004_Propofol\Modelling\head_models\sa_nyhead.mat')

neurons = {'L23E_oi24rpy1.swc';'L23I_oi38lbc1.swc';'L23I_oi38lbc1.swc';'L4E_53rpy1.swc'; ...
            'L4E_j7_L4stellate.swc';'L4E_j7_L4stellate.swc';'L4I_oi26rbc1.swc';'L4I_oi26rbc1.swc';
            'L5E_oi15rpy4.swc';'L5E_j4a.swc';'L5I_oi15rbc1.swc';'L5I_oi15rbc1.swc';'L6E_51-2a.CNG.swc'; ...
            'L6E_oi15rpy4.swc';'L6I_oi15rbc1.swc';'L6I_oi15rbc1.swc'};

neuron = neurons{5};
n=1;
propofol = 0;
[err,prints] = system(['python "E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\code\simulation_update20220519.py" ' int2str(propofol) ' ' neuron ' ' int2str(n)]);
outputs = split(prints,char(10));
fileIdcs = find(cellfun(@(x)~isempty(strfind(x,'E:')),outputs));

saveFile = outputs{fileIdcs(1)};
data = load(saveFile);
eeg1 = getEEG(data.Q',sa);
[p_baseline,f] = pmtm(eeg1,2,[],1e3*16);


t = []; g = []; s = []; syns= [];
for i = 1:length(data.synapses.times)
    n = length(data.synapses.times{i});
    g = [g;i*ones(n,1)];
    s = [s;(data.synapses.ei(i)=='i')*ones(n,1)];
    t0 = double(data.synapses.times{i}');
    if(any(and(t0>1319,t0<1321)))
        syns(end+1) = i;
    end
    t = [t;t0];
end

fig = figureNB;
    scatter(t,g,5,s,'filled')
    colormap([1,0,0;0,0,1])
    yticks([])
    gcaformat
    e0 = sum(data.synapses.ei=='e');
    line(get(gca,'xlim'),(e0+1.5)*[1,1],'color','k')
    xlim([1000,2000])
figureNB;
    axes('Position',[0,0,1,1]);
    x0 = data.morphology.coordinates(1,3:5);
    % x0 = mean(x0(x0(:,2)==1,3:5));
    for i = 1:length(data.morphology.segments)
        x = data.morphology.coordinates(data.morphology.segments{i},3:5)-x0;
        d = mean(data.morphology.coordinates(data.morphology.segments{i},6));
        plot3(x(:,1),x(:,2),x(:,3),'k','LineWidth',d); hold on;
    end
    set(gca,'DataAspectRatio',[1,1,1])
    axis off;
    % offset = 0*[50.6,2.366,-33.4];

    x1 = data.synapses.locations(data.synapses.ei=='e',:);%-[5,0,0];
    x2 = data.synapses.locations(data.synapses.ei=='i',:);%+[5,0,0];
    plot3(x1(:,1),x1(:,2),x1(:,3),'.','color','r','MarkerSize',4);
    plot3(x2(:,1),x2(:,2),x2(:,3),'.','color','b','MarkerSize',4);

    x = data.morphology.Vx; y = data.morphology.Vy; z = data.morphology.Vz;
    line(x(386,:),y(386,:),z(386,:),'LineWidth',5)

figureNB;
axes('Position',[0.15,0.1,0.79,0.35])
    y = data.V_soma; y(data.t<1000) = nan;
    plot(data.t,y,'color','k');
    xlim([925,2000]);
    line([950,980],[-65,-65],'LineWidth',1,'color','k')
    text(800,-65.3,{'V_{m} soma','-65 mV'},'FontSize',7,'HorizontalAlignment','left','color','k','VerticalAlignment','bottom')
    axis off;
    ylim([-70,-53]);
axes('Position',[0.15,0.55,0.79,0.35])
    % i = randi(size(data.V,1))
    y = data.V(386,:); y(data.t<1000) = nan;
    plot(data.t,y,'color','k');
    xlim([925,2000]);
    line([950,980],[-65,-65],'LineWidth',1,'color','k')
    text(800,-65.3,{'V_{m} dend','-65 mV'},'FontSize',7,'HorizontalAlignment','left','color','k','VerticalAlignment','bottom')
    axis off;
    ylim([-70,-53]);


figureNB;
axes('Position',[0.15,0.5,0.79,0.25])
    y = data.Q(1,:); y(data.t<1000) = nan;
    plot(data.t,y,'k');
    xlim([925,2000]);
    ylim([-10,10]);
    line([950,980],[0,0],'LineWidth',1,'color','k');
    text(800,-1.2,{'Px',['0 nA' char(956) 'm']},'FontSize',7,'HorizontalAlignment','left','color','k','VerticalAlignment','bottom')
    % text(940,0,['0 nA' char(956) 'm'],'FontSize',7,'HorizontalAlignment','right')
    % text(940,10,'Px','FontSize',7,'HorizontalAlignment','right')
    axis off;
axes('Position',[0.15,0.25,0.79,0.25])
    y = data.Q(2,:); y(data.t<1000) = nan;
    plot(data.t,y,'k');
    xlim([925,2000]);
    ylim([-10,10]);
    line([950,980],[0,0],'LineWidth',1,'color','k');
    text(800,-1.2,{'Py',['0 nA' char(956) 'm']},'FontSize',7,'HorizontalAlignment','left','color','k','VerticalAlignment','bottom')
    % text(940,0,['0 nA' char(956) 'm'],'FontSize',7,'HorizontalAlignment','right')
    % text(940,10,'Py','FontSize',7,'HorizontalAlignment','right')
    axis off;
axes('Position',[0.15,0,0.79,0.25])
    y = data.Q(3,:); y(data.t<1000) = nan;
    plot(data.t,y,'k'); 
    xlim([925,2000]);
    ylim([-10,10]);
    line([950,980],[0,0],'LineWidth',1,'color','k');
    text(800,-1.2,{'Pz',['0 nA' char(956) 'm']},'FontSize',7,'HorizontalAlignment','left','color','k','VerticalAlignment','bottom')
    % text(940,0,['0 nA' char(956) 'm'],'FontSize',7,'HorizontalAlignment','right')
    % text(940,10,'Pz','FontSize',7,'HorizontalAlignment','right')
    axis off;
    y0 = min(get(gca,'ylim'))*0.75;
    line([980,1100],[y0,y0],'color','k','LineWidth',1)
    text(1050,y0,'100 ms','FontSIze',7,'color','k','VerticalAlignment','top','HorizontalAlignment','center');




load('E:\Research_Projects\004_Propofol\Modelling\head_models\sa_nyhead.mat')
X = struct();
X.vertices = sa.cortex75K.vc;
X.faces= sa.cortex75K.tri;
figure;
    plot_mesh_brain(X);
    view([122,15]);
    fix_lighting;