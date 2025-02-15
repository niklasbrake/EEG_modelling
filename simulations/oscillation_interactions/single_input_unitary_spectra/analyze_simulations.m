function analyze_simulations(i)
addpath('/lustre04/scratch/nbrake/code/simulation_code');
folder = '/lustre04/scratch/nbrake/data/simulations/trend_peak_interaction2';

[sa,X] = network_simulation_beluga.getHeadModel;

locations = sa.cortex2K.in_from_cortex75K;

for i = 1:5
    for j = 1:5
        saveFolder = fullfile(folder,['rhythm' int2str(i-1)],['combo' int2str(j)]);
        load(fullfile(saveFolder,'simulation','simulation_data.mat'));
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
        f = [0.1:0.1:1e3]

        save(fullfile(saveFolder,'unitary_spectrum.mat'),'f','P');
    end
end