function changeEI_with_oscillation(i)
% baseFolder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\EI_oscillation';

baseFolder = '/lustre04/scratch/nbrake/data/simulations/EI_oscillation';
addpath('/lustre04/scratch/nbrake/code/simulation_code');

alpha = 2.^(-4:0);
% Parallelize with SLURM
% for i = 1:length(alpha)
for j = 1:5
    folder = fullfile(baseFolder,num2str(alpha(i)),['run' int2str(j)]);
    runSimulation(folder,alpha(i));
end
% end


function runSimulation(folder,alpha)
    % Initialize network
    network = network_simulation_beluga(folder);

    % Initialize post network
    nPostNeurons = 1;
    network = network.initialize_postsynaptic_network(nPostNeurons,1);

    % Presyanptic network parameters
    nPreNeurons = 30000;
    network.tmax = 10e3; % 2 seconds
    network.branchNo = 0;
    oscillation_frequency = 2;


    % eFiringRate = network.eFiringRate;
    % iFiringRate = network.iFiringRate;
    eFiringRate = 1;
    iFiringRate = 5;

    % Maintain mean firing rate
    b = @(x) network.eiFraction*(x-1)*eFiringRate/(1-network.eiFraction);

    network.eFiringRate = eFiringRate*alpha;
    network.iFiringRate = iFiringRate-b(alpha);

    simulatespikes_osc(nPreNeurons,network.tmax+100,network,oscillation_frequency);

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
