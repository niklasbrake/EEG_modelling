function aperiodic_sensitivity(iChunk)

% baseFolder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\aperiodic_interactions';
baseFolder = '/lustre04/scratch/nbrake/data/simulations/aperiodic_sensitivity';
addpath('/lustre04/scratch/nbrake/code/simulation_code');

fid = fopen(fullfile('/lustre04/scratch/nbrake/code/simulation_code','mod_files','default_parameters.json'));
str = char(fread(fid,inf)');
fclose(fid);
pars = jsondecode(str);

parameter_names = {'tauI',,'eL','lambdaE','m'};

P1 = [pars.iSynParams.tau2; ...
        pars.biophys_pars.pas_mem_pars.erev_leak; ...
        pars.eCellParams.firingRate; ...
        0];
P2 = [20;-45;2;0.98];

P = zeros(11,4);
P(1,:) = P1;
count = 2;
for i = 1:4
    for j = i:4
        P(count,:) = P1;
        P(count,i) = P2(i);
        P(count,j) = P2(j);
        count = count+1;
    end
end

M = 50;
V = zeros(32001,size(P,1)*M);
dipoles = zeros(32001,3,size(P,1)*M);
parameters = zeros(4,size(P,1)*M);
for i = 1:size(P,1);
    pars.iSynParams.tau2 = P(i,1);
    pars.biophys_pars.pas_mem_pars.erev_leak = P(i,2);
    pars.eCellParams.firingRate = P(i,3);
    mValue = P(i,4);

    for k = 1:M
        sampleID = M*(i-1) + k;
        folder = fullfile(baseFolder,'simulations',['sample' int2str(sampleID)]);
        dataFile = fullfile(folder,'simulation','simulation_data.mat');
        parameters(:,sampleID) = P(i,:);
        if(exist(dataFile))
            data = load(dataFile);
            V(:,sampleID) = data.V;
            dipoles(:,:,sampleID) = data.dipoles;
            time = data.time;
        else
            [time,V(:,sampleID),dipoles(:,:,sampleID)] = runSimulation(folder,pars,mValue);
        end
    end
end

save(fullfile(baseFolder,['analyzed_results_' int2str(iChunk) '.mat']),'time','V','dipoles','parameters','parameter_names');

end
function [time,V,dipoles] = runSimulation(folder,pars,mValue)
    % Initialize network
    network = network_simulation_beluga(folder,pars);

    % Initialize post network
    network = network.initialize_postsynaptic_network(1);

    % Presyanptic network parameters
    network.tmax = 2e3; % 2 seconds

    if(mValue==0)
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
    else
        nPreNeurons = 30e3;
        network.branchNo = mValue;
        network.simulatespikes_critplane(nPreNeurons,network.tmax);
        copyfile(fullfile(network.preNetwork,'spikeTimes.csv'),fullfile(network.preNetwork,'spikeTimesLong.csv'));
        x = csvread(fullfile(network.preNetwork,'locations.csv'));
        x1 = 2*pi*x(:,1);
        y1 = asin(2*x(:,2)-1)-pi/2;
        data = [(1:nPreNeurons)'-1,x1(:),y1(:)];
        dlmwrite(fullfile(network.preNetwork,'UMAP_embedding.csv'),data,'precision','%.4f');

        % Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
        network.form_connections(1);
    end

    % Simulate dipoles
    network = network.simulate();

    % Cleanup
    load(fullfile(network.savePath,'simulation_data.mat'));
end