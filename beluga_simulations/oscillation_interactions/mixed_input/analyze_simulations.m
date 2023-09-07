function analyze_simulations(i)
addpath('/lustre04/scratch/nbrake/code/simulation_code');
folder = '/lustre04/scratch/nbrake/data/simulations/mixed_input';

[sa,X] = network_simulation_beluga.getHeadModel;

% compute_spectrum(sa,fullfile(folder,'baseline'))
% compute_spectrum(sa,fullfile(folder,'high_aperiodic'))
% compute_spectrum(sa,fullfile(folder,'high_both'))
% compute_spectrum(sa,fullfile(folder,'high_oscillation'))
compute_spectrum(sa,fullfile(folder,'oscillation_only'))
compute_spectrum(sa,fullfile(folder,'aperiodic_only'))


end
function compute_spectrum(sa,folder)
    locations = sa.cortex2K.in_from_cortex75K;
    load(fullfile(folder,'simulation','simulation_data.mat'));
    m = size(dipoles,3);
    dp = [];
    for k = 1:m
        dp(:,:,k) = detrend(resample(dipoles(:,:,k),1e3,16e3),'constant');
    end
    eeg = zeros(size(dp,1),m*length(locations));
    G = zeros(1,m*length(locations));
    for k = 1:length(locations)
        eeg(:,m*(k-1)+1:m*k) = network_simulation_beluga.getEEG(dp,sa,locations(k));
        G(m*(k-1)+1:m*k) = 1:m;
    end
    psd = mypmtm(eeg,1e3,10);
    P = splitapply(@(x) mean(x,2),psd,G);
    f = [0.1:0.1:1e3];

    save(fullfile(folder,'unitary_spectrum.mat'),'f','P');
end