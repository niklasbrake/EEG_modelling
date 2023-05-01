function aperiodic_sensitivity(i)

% baseFolder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\aperiodic_interactions';
baseFolder = '/lustre04/scratch/nbrake/data/simulations/aperiodic_interactions';
addpath('/lustre04/scratch/nbrake/code/simulation_code');


EI0 = 1;
EI_delta = 2.^linspace(-2,0,5);

m0 = 0;
mValues = [0,0.63,0.86,0.95,0.98,0.99];

tau0 = 10;
tauValues = [10,15,20,25,30];

params = [EI_delta,EI0*ones(1,length(mValues)+length(tauValues)); ...
        m0*ones(1,length(EI_delta)),mValues,m0*ones(1,length(tauValues)); ...
        tau0*ones(1,length(EI_delta)+length(mValues)),tauValues];


% Parallelize with SLURM
for j = 1:20
    a = params(1,i);
    b = params(2,i);
    c = params(3,i);
    folder = fullfile(baseFolder,['combo' int2str(j)],['run' int2str(j)]);
    runSimulation(folder,a,b,c);
end

end
function runSimulation(folder,EI_delta,mValue,tau)
    % Initialize network
    network = network_simulation_beluga(folder);

    % Changes to EI balance (EI_delta = 2.^(-4:0);) while maintaing avg. firing rate
    eFiringRate_default = network.parameters.eCellParams.firingRate;
    iFiringRate_default = network.parameters.iCellParams.firingRate;
    b = @(x) network.parameters.eiFraction*(x-1)*eFiringRate_default/(1-network.parameters.eiFraction);
    network.parameters.eCellParams.firingRate = eFiringRate_default*EI_delta;
    network.parameters.iCellParams.firingRate = iFiringRate_default-b(EI_delta);

    % Changes to network crticality mValue = [0,0.63,0.86,0.95,0.98,0.99];
    network.branchNo = mValue;

    % Changes to syantpic time constant mValue = [0,0.63,0.86,0.95,0.98,0.99];
    network.parameters.iSynParams.tau2 = tau;

    % Initialize post network
    nPostNeurons = 1;
    network = network.initialize_postsynaptic_network(nPostNeurons,1);

    % Presyanptic network parameters
    nPreNeurons = 30000;
    network.tmax = 10e3; % 2 seconds
    oscillation_frequency = 2;

    simulatespikes_osc(nPreNeurons,network.tmax,network,oscillation_frequency);

    % Explicit mapping of plane onto sphere to save time (bypass UMAP; ignore seam)
    copyfile(fullfile(network.preNetwork,'spikeTimes.csv'),fullfile(network.preNetwork,'spikeTimesLong.csv'));
    x = csvread(fullfile(network.preNetwork,'locations.csv'));
    x1 = 2*pi*x(:,1);
    y1 = asin(2*x(:,2)-1)-pi/2;
    data = [(1:nPreNeurons)'-1,x1(:),y1(:)];
    dlmwrite(fullfile(network.preNetwork,'UMAP_embedding.csv'),data,'precision','%.4f');

    % Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
    network.form_connections(1);

    % Simulate dipoles
    network = network.simulate();
end