function aperiodic_sensitivity(i)

% baseFolder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\aperiodic_interactions';
baseFolder = '/lustre04/scratch/nbrake/data/simulations/passive_sensitivity';
addpath('/lustre04/scratch/nbrake/code/simulation_code');

tauIRange = [5:5:30];
tauERange = [1:0.5:3.5];
erevRange = [-70:5:-45];
eFiringRange = round(10.^linspace(-1,1,6),2,'significant');
iFiringRange = round(10.^linspace(0,log10(20),6),2,'significant');
g_leak = round(10.^linspace(-4.5,-2.5,6),2,'significant');

N(1) = length(tauIRange);
N(2) = length(tauERange);
N(3) = length(erevRange);
N(4) = length(eFiringRange);
N(5) = length(iFiringRange);
N(6) = length(g_leak);

fid = fopen(fullfile('/lustre04/scratch/nbrake/code/simulation_code','mod_files','default_parameters.json'));
str = char(fread(fid,inf)');
fclose(fid);
pars = jsondecode(str);

parameter = interp1(cumsum(N),1:length(N),i,'next','extrap');
iPar = i-sum(N(1:parameter-1))

switch parameter
case 1
    pars.iSynParams.tau2 = tauIRange(iPar);
    folder = fullfile(baseFolder,'tauI',num2str(pars.iSynParams.tau2));
case 2
    pars.eSynParams.tau2 = tauERange(iPar);
    folder = fullfile(baseFolder,'tauE',num2str(pars.eSynParams.tau2));
case 3
    pars.biophys_pars.pas_mem_pars.erev_leak = erevRange(iPar);
    folder = fullfile(baseFolder,'erev',num2str(pars.biophys_pars.pas_mem_pars.erev_leak));
case 4
    pars.eCellParams.firingRate = eFiringRange(iPar);
    folder = fullfile(baseFolder,'eFiring',num2str(pars.eCellParams.firingRate));
case 5
    pars.iCellParams.firingRate = iFiringRange(iPar);
    folder = fullfile(baseFolder,'iFiring',num2str(pars.iCellParams.firingRate));
case 6
    pars.biophys_pars.pas_mem_pars.g_leak = g_leak(iPar);
    folder = fullfile(baseFolder,'gleak',num2str(pars.biophys_pars.pas_mem_pars.g_leak));
end

runSimulation(folder,pars);

end
function runSimulation(folder,pars)
    % Initialize network
    network = network_simulation_beluga(folder,pars);

    % Initialize post network
    nPostNeurons = 100;
    ratios = [27,5,9,18,5,2,4,3,12,9,6];
    nrnType = [];
    for j = 1:length(ratios)
        nrnType = [nrnType;j*ones(ratios(j),1)];
    end
    network = network.initialize_postsynaptic_network(nPostNeurons,nrnType);

    % Presyanptic network parameters
    network.tmax = 100; % 2 seconds
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
end