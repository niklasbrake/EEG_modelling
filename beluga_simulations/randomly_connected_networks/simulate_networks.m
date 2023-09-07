folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\randomly_connected_networks';
% Initialize network
network = network_simulation_beluga(folder);

% Initialize post network
nPostNeurons = 2;
mTypes = 1;
network = network.initialize_postsynaptic_network(nPostNeurons,[1,1]);

% Presyanptic network parameters
N = 30e3;
network.tmax = 5e3; % 2 seconds
network.branchNo = 0.98;

% Randomly connected E-I network
[ids,ts,ei] = EI_crit(N,network.tmax);
csvwrite(fullfile(network.preNetwork,'multisynapse_IDs.csv'),-ones(N,1));
network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,network.spikingFile)

network = network.compute_presynaptic_correlations(network.spikingFile);
network.embed_presyanptic_neurons;

% Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
network = network.form_connections(1);

% Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
network.simulate();