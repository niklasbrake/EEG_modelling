folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw';
% Initialize network
network = network_simulation_beluga(fullfile(folder,'osc_2ndOrder','run1'));

% Initialize post network
nPostNeurons = 1;
network = network.initialize_postsynaptic_network(nPostNeurons,1);

% Presyanptic network parameters
nPreNeurons = 30000;
network.tmax = 10e3; % 2 seconds
network.branchNo = 0;
oscillation_frequency = 2;
% simulatespikes_osc(nPreNeurons,network.tmax+100,network,oscillation_frequency);
simulatespikes_2ndOrder(nPreNeurons,network,oscillation_frequency);

% Explicit mapping of plane onto sphere to save time (bypass UMAP; ignore seam)
copyfile(fullfile(network.preNetwork,'spikeTimes.csv'),fullfile(network.preNetwork,'spikeTimesLong.csv'));
x = csvread(fullfile(network.preNetwork,'locations.csv'));
x1 = 2*pi*x(:,1);
y1 = asin(2*x(:,2)-1)-pi/2;
data = [(1:nPreNeurons)'-1,x1(:),y1(:)];
dlmwrite(fullfile(network.preNetwork,'UMAP_embedding.csv'),data,'precision','%.4f');

% Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
network.form_connections(1);

network.parameters.iSynParams.tau2 = 30;
network.savePath = 'E:\Research_Projects\004_Propofol\data\simulations\raw\osc_2ndOrder\simulation_prop';
% Simulate dipoles
network = network.simulate();
