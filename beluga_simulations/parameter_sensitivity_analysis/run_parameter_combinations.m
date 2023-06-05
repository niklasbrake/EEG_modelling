function aperiodic_sensitivity(iChunk)

% baseFolder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\aperiodic_interactions';
baseFolder = '/lustre04/scratch/nbrake/data/simulations/parameter_sensitivity_analysis';
addpath('/lustre04/scratch/nbrake/code/simulation_code');

load(fullfile(baseFolder,'parameters.mat'));

fid = fopen(fullfile('/lustre04/scratch/nbrake/code/simulation_code','mod_files','default_parameters.json'));
str = char(fread(fid,inf)');
fclose(fid);
pars = jsondecode(str);

for i = 1:6
    strLHS{i} = ['pars.' sampled_parameters.Properties.VariableNames{i}];
end

chunk_size = 500;
N = size(sampled_parameters,1);
m = floor(N/chunk_size);
chunk = reshape(1:(m*chunk_size),chunk_size,m);

V = zeros(32001,chunk_size);
dipoles = zeros(32001,3,chunk_size);
for k = 1:chunk_size
    k2 = chunk(k,iChunk);

    pars = jsondecode(str);
    for i = 1:6
        strRHS = num2str(sampled_parameters{k2,i});
        eval([strLHS{i} '=' strRHS]);
    end
    mType = sampled_parameters{k2,7};

    folder = fullfile(baseFolder,'simulations',['sample' int2str(k2)]);
    dataFile = fullfile(folder,'simulation','simulation_data.mat');
    if(exist(dataFile))
        data = load(dataFile);
        V(:,k) = data.V;
        dipoles(:,:,k) = data.dipoles;
        time = data.time;
    else
        [time,V(:,k),dipoles(:,:,k)] = runSimulation(folder,pars,mType);
    end
    parameters(:,k) = sampled_parameters{k2,:}';
end

sample_idcs = chunk(:,iChunk);
save(fullfile(baseFolder,['analyzed_results_' int2str(iChunk) '.mat']),'time','V','dipoles','parameters','sample_idcs');

end
function [time,V,dipoles] = runSimulation(folder,pars,mType)
    % Initialize network
    network = network_simulation_beluga(folder,pars);

    % Initialize post network
    network = network.initialize_postsynaptic_network(1,mType);

    % Presyanptic network parameters
    network.tmax = 2e3; % 2 seconds
    N = network.getsynapsecount;
    nE = floor(network.getsynapsecount*network.parameters.eiFraction);
    nI = N-nE;

    nE_spikes = floor(nE*network.parameters.eCellParams.firingRate*network.tmax/1e3);
    nI_spikes = floor(nI*network.parameters.iCellParams.firingRate*network.tmax/1e3);

    ei = [zeros(nE,1);ones(nI,1)];
    ts = network.tmax*rand(nE_spikes+nI_spikes,1);
    ids = [randi(nE,nE_spikes,1);nE+randi(nI,nI_spikes,1)];
    network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,network.spikingFile);

    % Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
    network.form_connections(0);

    % Simulate dipoles
    network = network.simulate();

    % Cleanup
    load(fullfile(network.savePath,'simulation_data.mat'));
end