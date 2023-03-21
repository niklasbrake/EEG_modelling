function initialize_network(masterPath,branchNo)
% masterPath = '/lustre04/scratch/nbrake/data/simulations/raw/correlation_test';
addpath('/home/nbrake/aperiodic_EEG_modelling/simulations/functions');

tmax0 = 20e3;
tmax = 2e3;
nrnCount = 4;

network = network_simulation_beluga(masterPath);
network = network.initialize_postsynaptic_network(nrnCount);
network.tmax = tmax;
network.branchNo = branchNo;
network.simulatespikes_critplane(tmax0*1e-3);

copyfile(fullfile(network.preNetwork,'spikeTimes.csv'),fullfile(network.preNetwork,'spikeTimesLong.csv'));

x = csvread(fullfile(network.preNetwork,'locations.csv'));
x1 = 2*pi*x(:,1);
y1 = asin(2*x(:,2)-1);
data = uint32([(1:network.getsynapsecount)',x1(:),y1(:)]);
dlmwrite(fullfile(network.preNetwork,'UMAP_embedding.csv'),data,'precision',ceil(log10(N)));

network.save();