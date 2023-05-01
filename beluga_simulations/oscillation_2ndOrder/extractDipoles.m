function extractDipoles

baseFolder = '/lustre04/scratch/nbrake/data/simulations/osc_2ndOrder';

for i = 1:10
    folder = fullfile(baseFolder,['run' int2str(i)]);
    load(fullfile(folder,'simulation','simulation_data.mat'));
    base.time(:,i) = time;
    base.dipoles(:,i) = dipoles(:,3);

    load(fullfile(folder,'simulation_prop','simulation_data.mat'));
    prop.time(:,i) = time;
    prop.dipoles(:,i) = dipoles(:,3);
end

save(fullfile(baseFolder,'dipoles.mat'),'prop','base');