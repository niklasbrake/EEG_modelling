function combine_unitary_spectra

addpath('/lustre04/scratch/nbrake/code/simulation_code');
folder = '/lustre04/scratch/nbrake/data/simulations/oscillation_sensitivity';

P2 = zeros(size(P,1),size(P,2),20);
count = 1;
rhythm = [];
combo = [];
for i = 1:4
    for j = 1:5
        saveFolder = fullfile(folder,['rhythm' int2str(i-1)],['combo' int2str(j)]);
        load(fullfile(saveFolder,'unitary_spectrum.mat'),'f','P');
        P2(:,:,count) = P;
        count = count+1;
        rhythm(count) = i-1;
        combo(count) = j;
    end
end

saveFile = '/lustre04/scratch/nbrake/data/simulation_analyzed/oscillation_interactions.mat';
save(saveFile,'f','P2','rhythm','combo');
