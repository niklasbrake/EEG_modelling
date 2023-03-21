function runSimulations(folder)
% N = 1;
% network = network_simulation_beluga(char(folder));
% network = network.initialize_postsynaptic_network(N);
% network = network.initialize_presynaptic_network(0,2e3);
% network.form_connections;
% system(['python functions/prep_simulations.py ' network.postNetwork]);

addpath('/home/nbrake/aperiodic_EEG_modelling/simulations/functions');
load(fullfile(folder,'data.mat'));
correlationFile = network.getCorrelationFile;
spikingFileLong = fullfile(fileparts(correlationFile),'spikeTimesLong.csv');
spikingFile = network.spikingFile;
tmax = network.tmax;
tmax0 = 100e3;

% Set up new network
network = network_simulation_beluga(folder);
network.tmax = tmax;
network = network.initialize_postsynaptic_network(2);
network = network.setCorrelationFile(correlationFile);
network = network.setSavePath(fullfile(network.outputPath,'LFPy'));

% Sample from long spike train
[ids,ts,ei] = network_simulation_beluga.getprenetwork(spikingFileLong);
tmax0 = max(ts);
tmax = tmax+200;
t0 = rand*(tmax0-tmax);
idcs = find(and(ts>=t0,ts<t0+tmax));
network_simulation_beluga.save_presynaptic_network(ids(idcs),ts(idcs)-t0,ei,32770,spikingFile)
network = network.setSpikingFile(spikingFile);

% Simulate
network.form_connections(1);
network = network.simulate();
network.save();