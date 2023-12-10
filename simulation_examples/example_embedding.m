folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\example_embedding';

i = 55;
rep = 41;

load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis2\parameters.mat')

fid = fopen(fullfile('C:\Users\brake\Documents\GitHub\EEG_modelling\model','mod_files','default_parameters.json'));
str = char(fread(fid,inf)');
fclose(fid);

for j0 = 1:8
    strLHS{j0} = ['pars.' sampled_parameters.Properties.VariableNames{j0}];
end

% New parameters for each rep
k0 = 100*(i-1)+rep;
pars = jsondecode(str);
for j0 = 1:8
    strRHS = num2str(sampled_parameters{k0,j0});
    eval([strLHS{j0} '=' strRHS ';']);
end

fid = fopen('C:\Users\brake\Documents\GitHub\EEG_modelling\simulations\critical_dipole_correlation\dyad_morphologies.txt');
txt = textscan(fid,'%s%s','Delimiter',',');
mTypes = cat(2,txt{:,1},txt{:,2});

M = {'L23E_oi24rpy1','L23I_oi38lbc1','L4E_53rpy1','L4E_j7_L4stellate','L4I_oi26rbc1','L5E_j4a','L5E_oi15rpy4','L5I_oi15rbc1','L6E_51_2a_CNG','L6E_oi15rpy4','L6I_oi15rbc1'};
m_ID(1) = find(strcmp(M,mTypes{i,1}));
m_ID(2) = find(strcmp(M,mTypes{i,2}));

% Initialize network
network = network_simulation_beluga(folder);

% Initialize post network
network = network.initialize_postsynaptic_network(2,m_ID);

network.spikingFile = fullfile(network.preNetwork,'spikeTimesLong.csv');
% Presyanptic network parameters
nPreNeurons = 30e3;
network.tmax = 40e3;
network.branchNo = 0.98;
network.simulatespikes_critplane(nPreNeurons,network.tmax);
network.save();

network = network.compute_presynaptic_correlations(network.spikingFile);
network.embed_presyanptic_neurons;

network.tmax = 2e3;
network.spikingFile = fullfile(network.preNetwork,['spikeTimes' int2str(rep) '.csv']);
[ids,ts,ei] = network.resimulate_critplane;
        network_simulation_beluga.save_presynaptic_network(ids,ts,ei,30e3,network.spikingFile)
network.parameters = pars;
network.simulate;
