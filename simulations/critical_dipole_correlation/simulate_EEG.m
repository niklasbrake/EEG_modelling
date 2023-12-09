function simulate_EEG(masterFolder,arrayID,S)
addpath('/lustre04/scratch/nbrake/code/simulation_code');

baseFolder = '/lustre04/scratch/nbrake/data/simulations/parameter_sensitivity_analysis_2';
load(fullfile(baseFolder,'parameters.mat'));

fid = fopen(fullfile('/lustre04/scratch/nbrake/code/simulation_code','mod_files','default_parameters.json'));
str = char(fread(fid,inf)');
fclose(fid);


for j0 = 1:8
    strLHS{j0} = ['pars.' sampled_parameters.Properties.VariableNames{j0}];
end

for j = 1:10
    i = 10*(arrayID-1)+j;
    folder = fullfile(masterFolder,['run_' num2str(i,'%.3d')]);
    for rep = 1:100

        % New parameters for each rep
        k0 = 100*(i-1)+rep;
        pars = jsondecode(str);
        for j0 = 1:8
            strRHS = num2str(sampled_parameters{k0,j0});
            eval([strLHS{j0} '=' strRHS ';']);
        end

        % savefile = fullfile(folder,['simulation' int2str(rep) '_S=0'],'simulation_data.mat');
        % if(~exist(savefile))
            run_simulation(folder,rep,pars,S)
        % end
    end
end
end

function run_simulation(folder,rep,pars,S)

    load(fullfile(folder,'model.mat'));
    if(rep==1)
        network = network.form_connections(S);
    end

    network.tmax = 2e3;

    network.spikingFile = fullfile(network.preNetwork,['spikeTimes' int2str(rep) '.csv']);
    if(~exist(network.spikingFile))
        [ids,ts,ei] = network.resimulate_critplane;
        network_simulation_beluga.save_presynaptic_network(ids,ts,ei,30e3,network.spikingFile)
    end

    network.parameters = pars;
    network.simulate(['simulation_S=' sprintf('%0.1f',S) '_rep' int2str(rep) '.mat'],'True');
end