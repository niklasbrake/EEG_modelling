
fig = figureNB(8,5);
load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\network_criticality_spectra_(s=1).mat')

axes('Position',[0.18,0.38,0.27,0.47]);
    idcs = and(f>0,f<=100);
    clrs = clrsPT.sequential(length(Cq)+3);
    clrs = clrs(5:end,:);
    for i = 1:8
        % h(i) = plotwitherror(f(idcs),squeeze(P(idcs,:,i)),'M','LineWidth',1,'color',clrs(i,:)); hold on;
        h(i) = plot(f(idcs),mean(squeeze(P(idcs,i,:)),2),'LineWidth',1,'color',clrs(i,:)); hold on;
    end
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;
    ylim([1e-17,1e-12]);
    text(0.035,6e-12,'A','FontWeight','normal','fontsize',12);
    text(0.035,6e-12,'Random syn. config.','FontWeight','bold','fontsize',7,'HorizontalAlignment','left');

% load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\network_criticality_spectra_(s=0).mat')
% load('E:\Research_Projects\004_Propofol\data\simulations\raw\dyad_network_criticality\power_spectra.mat')
load('C:\Users\brake\Documents\temp\UMAP\power_spectra.mat');

axes('Position',[0.65,0.38,0.27,0.47]);
    idcs = and(f>0,f<=100);
    clrs = clrsPT.sequential(length(Cq)+3);
    clrs = clrs(5:end,:);
    for i = 1:8
        % h(i) = plotwitherror(f(idcs),squeeze(P(idcs,:,i)),'M','LineWidth',1,'color',clrs(i,:)); hold on;
        h(i) = plot(f(idcs),mean(P(idcs,:,i),2),'LineWidth',1,'color',clrs(i,:)); hold on;
    end
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;
    ylim([1e-17,1e-12]);
    text(0.035,6e-12,'B','FontWeight','normal','fontsize',12);
    text(0.035,6e-12,'Optimal syn. config.','FontWeight','bold','fontsize',7,'HorizontalAlignment','left');

m = [0,0.5,0.75,0.9,0.95,0.98,0.99,0.999];
axes('Position',[0,0.04,1,0.12]);
    for i = 1:8
        x = floor((i-1)/4);
        a = (1-x);
        b = (i-4*x)-1.25;
        text(b+0.25,a,['m = ' num2str(m(i),3)],'FontSize',7,'VerticalAlignment','middle');
        line([b,b+0.2],[a,a]-0.05,'color',clrs(i,:),'LineWidth',1.5);
    end
    xlim([-0.5,4]);
    ylim([-0.5,1.5]);
    % text(1.7,2.4,'Spike propogation','fontsize',7,'HorizontalAlignment','center')
    % text(-0.5,-1,'Synapse configuration:','fontsize',7)
    axis off;
