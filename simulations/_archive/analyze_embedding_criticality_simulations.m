[sa,X] = network_simulation.getHeadModel;
folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\dyad_network_criticality';

idx = randi(size(X.vertices,1),100,1);

h = waitbar(0);
F = dir(folder);
F = F(3:end);
for i = 1:length(F)
    h = waitbar(i/length(F));
    folder2 = fullfile(folder,F(i).name);
    F2 = dir(folder2);
    F2 = F2(3:end);
    for k = 1:length(F2)
        load(fullfile(folder2,F2(k).name,'s=0','data.mat'));
        network = network.importResults;
        for k = 1:length(network.neurons)
            [~,pt1,pt2] = fileparts(network.neurons(k).sim);
            network.neurons(k).sim = fullfile(folder2,F2(k).name,'s=0','LFPy',[pt1 pt2]);
        end
        for j = 1:length(idx)
            eeg = network_simulation.getEEG(network.results.dipoles,sa,idx(j));
            idcs = find(network.results.t<network.results.t(end)-100);
            tempQ = corr(eeg(idcs,:));
            temp(j) =  mean(tempQ(tempQ~=1));
            [Ptemp(:,j),f] = pmtm(detrend(sum(eeg,2)),2,[],1e3*16);
        end
        Cq(i,k) = mean(temp);
        P(:,i,k) = mean(Ptemp,2);
    end
end

h = waitbar(0);
F = dir(folder);
F = F(3:end);
for i = 1:length(F)
    h = waitbar(i/length(F));
    folder2 = fullfile(folder,F(i).name);
    F2 = dir(folder2);
    F2 = F2(3:end);
    for k = 1:length(F2)
        load(fullfile(folder2,F2(k).name,'s=1','data.mat'));
        for l = 1:length(network.neurons)
            [~,pt1,pt2] = fileparts(network.neurons(l).sim);
            network.neurons(l).sim = fullfile(folder2,F2(k).name,'s=1','LFPy',[pt1 pt2]);
        end
        network = network.importResults;
        for j = 1:length(idx)
            eeg = network_simulation.getEEG(network.results.dipoles,sa,idx(j));
            idcs = find(network.results.t<network.results.t(end)-100);
            tempQ = corr(eeg(idcs,:));
            temp(j) =  mean(tempQ(tempQ~=1));
            % [Ptemp(:,j),f] = pmtm(detrend(sum(eeg,2)),2,[],1e3*16);
        end
        Cq_NULL(i,k) = mean(temp);
        % P(:,i,k) = mean(Ptemp,2);
    end
end




m = [0,0.5,0.75,0.9,0.95,0.98,0.99,0.999];
clrs = clrsPT.sequential(length(m)+4);
clrs = clrs(5:end,:);
figureNB(7.2,3.2);
subplot(1,2,1);
    plot(1./(1-m),mean(Cq,2),'color','k','LineWidth',1);
    hold on
    line([1,1000],[0,0],'color',[1,1,1]*0.7,'LineWidth',1)
    errorbar(1./(1-m),nanmean(Cq_NULL'),nanstd(Cq_NULL'),'.','Color',[0.7,0.7,0.7],'MarkerSize',10);
    for i = 1:length(m)
        h(i) = errorbar(1./(1-m(i)),mean(Cq(i,:)'),std(Cq(i,:)'),'.','Color',clrs(i,:),'MarkerSize',10);
    end
    set(gca,'xscale','log')
    gcaformat
    xlabel('1/(1-m)')
    ylabel('Dipole correlation')
    xticks([1,10,100,1000])
    xlim([1,1000])
    ylim([-0.05,0.8])
    yticks([0,0.4,0.8]);
subplot(1,2,2);
    for i = 1:length(m)
        plot(f,mean(P(:,i,:),3),'LineWidth',1,'color',clrs(i,:));
        hold on;
    end
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.25,150])
    xticks([1,10,100])
    xticklabels([1,10,100])
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    yticks([1e-18,1e-15,1e-12])
    ylim([1e-18,1e-12]);
    gcaformat


figureNB
for j = 1:length(m)/2
    idcs = (2*(j-1)+1:2*j)';
    ax = axes('Position',[0,0,1,1]);
    h = plot(nan,nan);
    h(1) = plot(nan,nan,'LineWidth',1,'color',clrs(idcs(1),:)); hold on;
    h(2) = plot(nan,nan,'LineWidth',1,'color',clrs(idcs(2),:));
    L = legend(h,strcat(repmat('m = ',[length(m(idcs)),1]),num2str(m(idcs)',3)));
    L.ItemTokenSize = [10,5];
    L.Box = 'off';
    L.Position(1) = (j-1)/5+0.02;
    ax.Color = 'none';
    gcaformat
end