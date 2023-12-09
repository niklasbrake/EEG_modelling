function combine_unitary_spectra

addpath('/lustre04/scratch/nbrake/code/simulation_code');
baseFolder = '/lustre04/scratch/nbrake/data/simulations/mixed_input';

[sa,X] = network_simulation_beluga.getHeadModel;


inputType = {'baseline', 'high_aperiodic', 'high_both', 'high_oscillation','oscillation_only','aperiodic_only'};

rhythm = [];
combo = [];
for i = 1:length(inputType)
    saveFolder = fullfile(baseFolder,inputType{i});
    load(fullfile(saveFolder,'simulation','simulation_data.mat'));
    m = size(dipoles,3);
    dp = [];
    for k = 1:m
        dp(:,:,k) = detrend(resample(dipoles(:,:,k),1e3,16e3),'constant');
    end

    load(fullfile(saveFolder,'unitary_spectrum.mat'),'f','P');
    if(i==1)
        P2 = zeros(size(P,1),size(P,2),4);
        DP = zeros(size(dp,1),size(dp,2),size(dp,3),4);
    end
    DP(:,:,:,i) = dp;
    P2(:,:,i) = P;
end


saveFile = '/lustre04/scratch/nbrake/data/simulation_analyzed/trend_peak_interaction_mixed.mat';
save(saveFile,'f','P2','DP','inputType');
