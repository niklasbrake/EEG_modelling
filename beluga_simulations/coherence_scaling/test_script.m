
[sa,X] = network_simulation_beluga.getHeadModel;
locations = randi(size(sa.cortex75K.vc,1),100,1); % Random location
load('E:\Research_Projects\004_Propofol\data\simulations\raw\peak_trend_sensitivity\trend_peak_interaction.mat')


figureNB(19,8);
idx = randi(size(sa.cortex75K.vc,1));
for i = 1:5
    eeg = network_simulation_beluga.getEEG(DP(:,:,:,5*(i-1)+1),sa,idx);
    % eeg = randn(10001,50);
    Y = eeg(500:end,:);
    P0 = 30*mean(pmtm(Y,2,[],1e3),2);
    PP = zeros(size(P0,1),200);
    for k = 1:size(PP,2)
        idcs = randperm(50,30);
        Y = eeg(500:end,idcs);
        [PP(:,k),f2] = pmtm(sum(Y,2),2,[],1e3);
        % PP1(:,k) = sum(pmtm(Y,2,[],1e3),2);
    end
    subplot(5,1,i);
    plot(f2,P0);
    hold on;
    plot(f2,mean(PP,2))
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    drawnow;
end

figureNB;
    plot(f2,(mean(PP,2)-P0)./P0/39);
    set(gca,'xscale','log')









m = 10;

P0 = zeros(2049,size(dipoles,4),m);
P1 = zeros(2049,size(dipoles,4),m);
for k = 1:size(dipoles,4)
    for i = 1:m
        waitbar(i/m)
        eeg = network_simulation_beluga.getEEG(dipoles(:,:,:,k),sa,randi(size(sa.cortex75K.vc,1)));
        P0(:,k,i) = sum(pmtm(eeg,2,[],1e3),2);
        [P1(:,k,i),f1] = pmtm(sum(eeg,2),2,[],1e3);
    end
end
figureNB;
    plot(f1,mean(mean(P0,3),2))
    hold on;
    plot(f1,mean(mean(P1,3),2))


figureNB;
    plot(f1,(mean(mean(P1,3),2)-mean(mean(P0,3),2))./mean(mean(P0,3),2))
    set(gca,'xscale','log')