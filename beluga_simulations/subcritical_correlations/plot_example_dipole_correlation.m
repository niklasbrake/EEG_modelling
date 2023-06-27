[sa,X] = network_simulation_beluga.getHeadModel;
folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\dyad_network_criticality\m=0.98';
idx = randi(size(sa.cortex75K.vc,1));
for k = 1:2
    X = csvread(fullfile(folder,['run0'],'\s=0\LFPy\',['cell0000' int2str(k) '.csv']),1,0);
    % eeg(:,k,j) = network_simulation_beluga.getEEG(X(2:end-16e2,2:4),sa,idx);
    dipoles(:,:,k) = X(2:end-16e2,2:4);
end

time = X(2:end-16e2,1);
eeg = network_simulation_beluga.getEEG(dipoles(:,:,:,1),sa,idx);

figureNB(5.25,1)
X = eeg(3e4:6e4,:);
axes('Position',[0,0,1,1])
    plot(X(:,1),'color','r')
    hold on;
    plot(X(:,2),'color','k')
    line([38e3-250*16,38e3],[-1.5,-1.5]*1e-6,'color','k','LineWidth',1)
    line([38e3-250*16,38e3-250*16],[-1.25,-0.25]*1e-6,'color','k','LineWidth',1)
    xlim([1,40e3]);
    ylim([-2e-6,0])
    axis off;