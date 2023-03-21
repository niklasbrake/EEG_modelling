function network = simulateNetwork(networkPath,spikingFile,savePath,tmax)

pyPath = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/code';
neuronSimulation = 'computeDipole.py';
computeDipole = fullfile(pyPath,neuronSimulation);

[mFiles,synapseFiles] = import_postsyanptic_network(networkPath);
saveOutput = fullfile(savePath,'output.mat');

mkdir(fullfile(savePath,'LFPy'));


[err,prints] = system(['python "' computeDipole '" ' params]);
if(err)
    error(prints);
end
outputFiles = split(prints,char(10));
idx = find(cellfun(@(x)~isempty(strfind(x,'E:')),outputFiles));
for j = 1:length(idx)
    outputFile = outputFiles{idx(j)};
    network(j) = simneuron(outputFile,neuron)
end
save(saveOutput,network);
