
tmax = 4; % seconds
propofol = 0;
n = 100; % Total size of simulation
similarityIndex = 0.5;
dynamics = 'critical';
% dynamics = 'asynchronous';

pyPath = 'E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/code';
neuronSimulation = 'simulation_synapse_similarity.py';
C = clock; today = sprintf('%02d',C(1:3));
savePath = ['E:/Research_Projects/004_Propofol/Modelling/neuron_simulations/data/simulations/compare_networks/' dynamics];
mkdir(savePath);
synapseFolder = fullfile(savePath,'network');
mkdir(synapseFolder);
simulationFolder = fullfile(savePath,'LFPy');
mkdir(simulationFolder);

neurons = {'L23E_oi24rpy1';'L23I_oi38lbc1';'L23I_oi38lbc1';'L4E_53rpy1'; ...
            'L4E_j7_L4stellate';'L4E_j7_L4stellate';'L4I_oi26rbc1';'L4I_oi26rbc1';
            'L5E_oi15rpy4';'L5E_j4a';'L5I_oi15rbc1';'L5I_oi15rbc1';'L6E_51_2a_CNG'; ...
            'L6E_oi15rpy4';'L6I_oi15rbc1';'L6I_oi15rbc1'};
noOccurances = [26.8,3.2,4.3,9.5,9.5,9.5,5.6,1.5,4.9,1.3,0.6,0.8,14,4.6,1.9,1.9]/100*n; % Hagen et al., 2016
expectedNoSynapses = [10717,5463,5463,8244,4968,4968,2713,2713,16384,6578,5904,5904,6516,7110,5904,5904]; % 1.15 syn/um
synapsesPerConnection = 4.5; % Markram et al. 2014
totalConnections = sum(expectedNoSynapses.*noOccurances)/synapsesPerConnection; % Total number of  connections in EEG model
% Average density of connections per presyanptic neuron per postsynaptic neuron (Markram et al., 2014)
rho = 254/31000;
i = 0; N = 0;
while N<totalConnections
    i = i+1;
    N = N+1+poissrnd(rho*n); % Random number of connections/presyanptic neuron
end
%%%%%%%%%%%%%%%%%%%%%%%%%% Define connections  %%%%%%%%%%%%%%%%%%%%%%%%
load('E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\cortical_column_Hagen\segment_areas.mat')

mTypes = unique(neurons);
mTypeCount = accumarray(findgroups(neurons),noOccurances);
occuranceCDF = cumsum(mTypeCount)/sum(mTypeCount);
randM = interp1(occuranceCDF,1:length(occuranceCDF),rand(n,1),'next','extrap');
M = histcounts(randM,0.5:1:length(mTypes)+0.5);
nSyns = splitapply(@mean,expectedNoSynapses(:),findgroups(neurons));

% Define geometric correlation among synapse connections
X = [];
for i = 1:length(mTypes)
    neuron = mTypes{i};
    x = nrnSegSA.(neuron).pos;
    X = [X;repmat(x,[ceil(mTypeCount(i)),1])];
end
X = X./vecnorm(X,2,2);

nCommon = ceil(similarityIndex*mean(expectedNoSynapses));
posCommon = X(randi(size(X,1),nCommon,1),:);
preCommon = randsample(N,nCommon);

specialPrePool = setdiff(1:N,preCommon);

nrns = struct();
for i = 1:length(mTypes)
    neuron = mTypes{i};
    reps = M(i);
    SA = nrnSegSA.(neuron).area;
    segPos = nrnSegSA.(neuron).pos;
    segPos =segPos./vecnorm(segPos,2,2);
    D = zeros(size(segPos,1),size(posCommon,1));
    for j = 1:size(posCommon,1)
        D(:,j) = sum(segPos.*posCommon(j,:),2);
    end
    [~,segCommon] = max(D);
    sa = cumsum(SA)/sum(SA);
    for k = 1:reps
        simID = [neuron '_N' int2str(k)];
        filename = fullfile(synapseFolder,[simID '.csv']);
        fid = fopen(filename,'w');
        nTotal = poissrnd(nSyns(i));
        preSpecial = randsample(specialPrePool,nTotal-nCommon);
        segSpecial = interp1(sa,1:length(sa),rand(nTotal-nCommon,1),'next','extrap');
        pre = [preSpecial(:);preCommon(:)];
        seg = [segSpecial(:);segCommon(:)];
        for j = 1:nTotal
            fprintf(fid,'%s,%s\n',int2str(seg(j)),int2str(pre(j)));
        end
        fclose(fid);
    end
end


%%%%%%%%%%%%%%%%%% Compute presyanptic network %%%%%%%%%%%%%%%%%%
if(strcmp(dynamics,'critical'))
    m = 0.98;
else
    m = 0;
end
spikingFile = fullfile(savePath,'network_spikes.csv');
[neuronID,spikeTime,ei,B] = simulatespikes(tmax,N,m);
fid = fopen(spikingFile,'w');
[uM,~,nTemp] = unique(neuronID(:));
I = accumarray(nTemp,1:numel(nTemp),[],@(x){sort(x)});
UC = 0;
flop = rand(max(neuronID),1);
for i = 1:max(neuronID)
    if(uM(i-UC)==i)
        ts = spikeTime(I{i-UC});
    else
        ts = [];
        UC = UC+1;
    end
    if(ei(i)==0)
        if(flop(i)<0.5)
            target = 'eA';
        elseif(flop(i)<=1)
            target = 'eB';
        end
    else
        if(flop(i) < 0.5)
            target = 'iA';
        elseif(flop(i) <= 0.9)
            target = 'iB';
        elseif(flop(i) <= 1)
            target = 'iS';
        end
    end
    s = sprintf('%0.1f,',sort(ts));
    fprintf(fid,'%s,%s\n',target,s(1:end-1));
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Simuate dipoles  %%%%%%%%%%%%%%%%%%%%%%%%
saveOutput = fullfile(savePath,'output.mat');
for i = 1:length(mTypes)
    neuron = mTypes{i};
    reps = M(i);
    if(reps==0)
        continue
    end
    [err,prints] = system(['python "' fullfile(pyPath,neuronSimulation) '" ' int2str(propofol) ' ' neuron ' ' int2str(reps) ' ' savePath ' ' int2str(tmax*1e3)]);
    outputFiles = split(prints,char(10));
    idx = find(cellfun(@(x)~isempty(strfind(x,'E:')),outputFiles));
    for j = 1:length(idx)
        outputFile = outputFiles{idx(j)};
        simID = [neuron '_N' int2str(j)];
        assignin('base',simID,simneuron(outputFile,neuron))
        if(exist(saveOutput))
            save(saveOutput,simID,'-append');
        else
            save(saveOutput,simID);
        end
    end
end


