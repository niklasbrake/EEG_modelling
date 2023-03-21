function runSimulations(network,repNo)

addpath('/home/nbrake/aperiodic_EEG_modelling/simulations/functions');

network.spikingFile = fullfile(network.preNetwork,['spikeTimes_' int2str(repNo) '.csv']);
network.savePath = [network.savePath '_' int2str(repNo)];

[ids,ts,ei] = network.getprenetwork(fullfile(fileparts(network.getCorrelationFile),'spikeTimesLong.csv'));
tmax = network.tmax;
tmax0 = max(ts);
t0 = rand*max(tmax0-tmax,0);
idcs = find(and(ts>=t0,ts<t0+tmax));
N = network.getsynapsecount;
network_simulation_beluga.save_presynaptic_network(ids(idcs),ts(idcs)-t0,ei,N,network.spikingFile)


network.form_connections(0.5);
network = network.simulate();