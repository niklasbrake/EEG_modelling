folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw';
% Initialize network
network = network_simulation_beluga(fullfile(folder,'EI_oscillations','WTising'));

% Initialize post network
nPostNeurons = 1;
network = network.initialize_postsynaptic_network(nPostNeurons,1);

% Presyanptic network parameters
N = 30000;
network.tmax = 2e3; % 2 seconds
network.branchNo = 0;

% Simulate spike trains
[ids,ts,ei,x] = simulatespikes_avaETosc(network,0.987,0.01,N);
network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,network.spikingFile)
csvwrite(fullfile(network.preNetwork,'multisynapse_IDs.csv'),-ones(N,1));

% Explicit mapping of plane onto sphere to save time (bypass UMAP; ignore seam)
copyfile(fullfile(network.preNetwork,'spikeTimes.csv'),fullfile(network.preNetwork,'spikeTimesLong.csv'));
x1 = 2*pi*x(:,1);
y1 = asin(2*x(:,2)-1)-pi/2;
data = [(1:N)'-1,x1(:),y1(:)];
dlmwrite(fullfile(network.preNetwork,'UMAP_embedding.csv'),data,'precision','%.4f');

% Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
network.form_connections(1);

% Simulate dipoles
network = network.simulate();