function combine_unitary_spectra

addpath('/lustre04/scratch/nbrake/code/simulation_code');
folder = '/lustre04/scratch/nbrake/data/simulations/trend_peak_interaction2';

count = 1;
rhythm = [];
combo = [];

for i = 1:5
    for j = 1:5
        saveFolder = fullfile(folder,['rhythm' int2str(i-1)],['combo' int2str(j)]);
        load(fullfile(saveFolder,'simulation','simulation_data.mat'));
        m = size(dipoles,3);
        dp = [];
        for k = 1:m
            dp(:,:,k) = detrend(resample(dipoles(:,:,k),1e3,16e3),'constant');
        end

        load(fullfile(saveFolder,'unitary_spectrum.mat'),'f','P');
        if(count==1)
            P2 = zeros(size(P,1),size(P,2),25);
            DP = zeros(size(dp,1),size(dp,2),size(dp,3),25);
        end
        DP(:,:,:,count) = dp;
        P2(:,:,count) = P;
        count = count+1;
        rhythm(count) = i-1;
        combo(count) = j;
    end
end

saveFile = '/lustre04/scratch/nbrake/data/simulation_analyzed/trend_peak_interaction.mat';
save(saveFile,'f','P2','rhythm','combo','DP');
