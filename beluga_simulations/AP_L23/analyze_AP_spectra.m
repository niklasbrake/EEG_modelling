passive = load('E:\Research_Projects\005_Aperiodic_EEG\data\simulations\active\LFPy_passive\simulation_data.mat')
active = load('E:\Research_Projects\005_Aperiodic_EEG\data\simulations\active\LFPy_active\simulation_data.mat')

idcs = find(passive.time>100);

passive.time = passive.time(idcs);
passive.V = passive.V(idcs,:);
passive.dipoles = passive.dipoles(idcs,:,:);

active.time = active.time(idcs);
active.V = active.V(idcs,:);
active.dipoles = active.dipoles(idcs,:,:);

[sa,X] = network_simulation_beluga.getHeadModel;
location = randi(size(sa.cortex75K.vc,1)); % Random location

passive.eeg = network_simulation_beluga.getEEG(passive.dipoles,sa,location);
[passive.psd,f] = pmtm(detrend(passive.eeg),2,[],16e3);

active.eeg = network_simulation_beluga.getEEG(active.dipoles,sa,location);
[active.psd,f] = pmtm(detrend(active.eeg),2,[],16e3);


figureNB;
    plot(f,mean(active.psd,2))
    hold on;
    plot(f,mean(passive.psd,2))
    set(gca,'xscale','log')
    set(gca,'yscale','log')