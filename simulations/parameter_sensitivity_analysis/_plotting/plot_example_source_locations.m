% Load NY head model
[sa,X] = network_simulation_beluga.getHeadModel;
idcs = sa.cortex2K.in_from_cortex75K;
czIdx = find(strcmp(sa.clab_electrodes,'Cz'));

% Note: different simulation from published figure
load(fullfile(dataFolder,'simulations','example_simulation','simulation','simulation_data.mat'));

% Compute location-averaged spectrum
N = size(sa.cortex75K.vc,1);
dp = resample(dipoles,250,16e3); % resample to 250 Hz
eeg = zeros(size(dp,1),N);
for i = 1:N
    eeg(:,i) = network_simulation_beluga.getEEG(dp,sa,i);
end
[psd,f] = pmtm(eeg,2,[],250);

iExample = [32918,37229];

x1 = [X.vertices(iExample(1),:),sa.cortex75K.normals(iExample(1),:)];
x2 = [X.vertices(iExample(2),:),sa.cortex75K.normals(iExample(2),:)];

blue = [0,1,1]/2;
red = [1,0,1]/1.5;

fig = figure('color','w','units','centimeters');
fig.Position(3:4) = [6.5,3.3];
axes('Position',[0.0,-0.2,0.46,1.54]);
    trisurf(X.faces, X.vertices(:,1), X.vertices(:,2), X.vertices(:,3),...
    'FaceLighting','gouraud','FaceVertexCData',0,'EdgeColor','none');
    hold on
    plot3(X.vertices(idcs,1),X.vertices(idcs,2),X.vertices(idcs,3),'.k','MarkerSize',2);
    plot3(X.vertices(iExample(1),1),X.vertices(iExample(1),2),X.vertices(iExample(1),3),'.','color',blue,'MarkerSize',15);
    plot3(X.vertices(iExample(2),1),X.vertices(iExample(2),2),X.vertices(iExample(2),3),'.','color',red,'MarkerSize',15);
    plot3(sa.locs_3D(czIdx,1),sa.locs_3D(czIdx,2),sa.locs_3D(czIdx,3),'.k','MarkerSize',10);
    view([122,15]);
    axis tight equal off
    camlight headlight
    material dull
    colormap(flip(gray(10)))
    caxis([0,1]);

    L = drawline_3Dplot(x1(1:3),pi/4,30)
    L.Color = blue;
    text(L.XData(2),L.YData(2),L.ZData(2),'Location A','fontsize',6,'HorizontalAlignment','left','VerticalAlignment','bottom','color',L.Color)

    L = drawline_3Dplot(x2(1:3),-3*pi/4,40);
    L.Color = red;
    text(L.XData(2),L.YData(2),L.ZData(2),'Location B','fontsize',6,'HorizontalAlignment','center','VerticalAlignment','top','color',L.Color)

    L = drawline_3Dplot(sa.locs_3D(czIdx,1:3),pi/2,20);
    pos = [L.XData(2),L.YData(2),L.ZData(2)];
    delete(L);
    text(pos(1),pos(2),pos(3),'Electrode','fontsize',6,'HorizontalAlignment','center','VerticalAlignment','middle')

    L = drawline_3Dplot([0,0,0],-pi/2,90);
    pos = [L.XData(2),L.YData(2),L.ZData(2)];
    delete(L);
    text(pos(1),pos(2),pos(3),'Neuron placement','fontsize',7,'HorizontalAlignment','center','VerticalAlignment','top')


axes('Position',[0.68,0.3,0.28,0.57]);
    plot(f,psd(:,iExample(1)),'color',blue,'LineWidth',0.5); hold on;
    plot(f,psd(:,iExample(2)),'color',red,'LineWidth',0.5); hold on;
    plot(f,mean(psd,2),'k','LineWidth',1)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,150]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xl = xlabel('Frequency (Hz)');
    xl.Position(2) = 1.6055e-18;
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    ylim([1e-17,1e-14])
    yticks([1e-17,1e-14])
    gcaformat;
    title('EEG spectrum','fontweight','normal','fontsize',7)
    text(0.75,2^3.4*1e-17,'Loc. A','FontSize',6,'color',blue)
    text(0.75,2^2.2*1e-17,'Loc. B','FontSize',6,'color',red)
    text(0.75,2*1e-17,'Cortex avg.','FontSize',6,'color','k')