[sa,X] = network_simulation_beluga.getHeadModel;
folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\synthetic_spikes';

locations = randi(size(sa.cortex75K.vc,1),100,1); % Random location

C = zeros(6,1);
for i = 1:6
    if(i==3)
        C(i)=nan;
        continue;
    end
    load(fullfile(folder,['mCombo13_R' int2str(i)],'simulation_spaceshuffle\simulation_data.mat'))
    for k = 1:100
        eeg = network_simulation_beluga.getEEG(dipoles,sa,locations(k));
        C(i) = C(i) + corr(eeg(:,1),eeg(:,2))/100;
    end
end


rRange = 10.^linspace(-4,-2,6);

figureNB(3.36,3.05);
    plot(rRange([1,2,4,5,6]),C2([1,2,4,5,6]),'.-k','MarkerSize',10);
    hold on;
    % plot(rRange([1,2,4,5,6]),C([1,2,4,5,6]),'.-','color',[0.6,0.6,0.6],'MarkerSize',10);
    xlim([-0.0005,0.01])
    ylim([-0.02,0.2])
    ylabel('\rho_{max}')
    set(gca,'xscale','log')
    xlabel('Spike correlation (R^2)')
    gcaformat