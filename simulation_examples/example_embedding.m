folder = fullfile(pwd,'example_embedding');

% Initialize network
network = network_simulation_beluga(folder);

% Initialize post network
network = network.initialize_postsynaptic_network(2,[1,2]);

network.spikingFile = fullfile(network.preNetwork,'spikeTimesLong.csv');
% Presyanptic network parameters
nPreNeurons = 30e3;
network.tmax = 40e3;
network.branchNo = 0.98;
network.simulatespikes_critplane(nPreNeurons,network.tmax);
network.save();

network = network.compute_presynaptic_correlations(network.spikingFile);
network.embed_presyanptic_neurons;

network.tmax = 4e3;
network.spikingFile = fullfile(network.preNetwork,'spikeTimes.csv');
[ids,ts,ei] = network.resimulate_critplane;
        network_simulation_beluga.save_presynaptic_network(ids,ts,ei,nPreNeurons,network.spikingFile)
network.parameters = pars;
network.simulate;
