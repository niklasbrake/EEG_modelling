% function runSimulations(masterPath)

% masterPath='/lustre04/scratch/nbrake/data/simulations/raw/correlation_test';

addpath('/home/nbrake/aperiodic_EEG_modelling/simulations/functions');
load(fullfile(masterPath,'data.mat'));
[ids,ts,ei] = network.getprenetwork(fullfile(fileparts(network.getCorrelationFile),'spikeTimesLong.csv'));
tmax = network.tmax;
tmax0 = max(ts);
t0 = rand*max(tmax0-tmax,0);
idcs = find(and(ts>=t0,ts<t0+tmax));
N = network.getsynapsecount;
network_simulation_beluga.save_presynaptic_network(ids(idcs),ts(idcs)-t0,ei,N,network.spikingFile)

network.form_connections(1);
network = network.simulate();