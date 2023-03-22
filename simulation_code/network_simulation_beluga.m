classdef network_simulation_beluga

    properties
        tmax
        branchNo
        neurons
        outputPath
        postNetwork
        preNetwork
        spikingFile
        savePath
    end

    properties (Constant)
        resourceFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';
        functionFolder = fileparts(mfilename('fullpath'));
        eiFraction = 0.85;
        eFiringRate = 0.5; % Hz
        iFiringRate = 2.5; % Hz
        eMulti = 1;
        iMulti = 1;
        % eMulti = 3.6; % Syanpses/connections Markram et al.
        % iMulti = 13.9; % Syanpses/connections Markram et al.
    end

    properties (Access = private)
        morphologyPath = fullfile(network_simulation_beluga.resourceFolder,'cortical_column_Hagen','swc');
        mTypeSegmentationData = fullfile(network_simulation_beluga.resourceFolder,'cortical_column_Hagen','morphology_segmentations.mat');
        simulateFunction = fullfile(network_simulation_beluga.functionFolder,'compute_network_dipoles_beluga.py');
        convertFilesFunction = fullfile(network_simulation_beluga.functionFolder,'np2mat.py');
        correlationFunction = fullfile(network_simulation_beluga.functionFolder,'compute_tiling_correlation.exe');
        embeddingFunction = fullfile(network_simulation_beluga.functionFolder,'embed_data.py');
        neuronTypes = {'L23E_oi24rpy1';'L23I_oi38lbc1';'L23I_oi38lbc1';'L4E_53rpy1';'L4E_j7_L4stellate';'L4E_j7_L4stellate';'L4I_oi26rbc1';'L4I_oi26rbc1';'L5E_oi15rpy4';'L5E_j4a';'L5I_oi15rbc1';'L5I_oi15rbc1';'L6E_51_2a_CNG';'L6E_oi15rpy4';'L6I_oi15rbc1';'L6I_oi15rbc1'};
        nrnAbundance = [26.8,3.2,4.3,9.5,9.5,9.5,5.6,1.5,4.9,1.3,0.6,0.8,14,4.6,1.9,1.9]/100
        dendriteLengths = [9319,4750,4750,7169,4320,4320,2359,2359,14247,5720,5134,5134,5666,6183,5134,5134];
        activeConductances = false
        inSynParamChanges = 0
        synapseCount
        morphoCounts
        correlationFile
    end

    methods

        function obj = network_simulation_beluga(outputPath)
            warning('off','MATLAB:MKDIR:DirectoryExists')
            if(nargin>0)
                obj.outputPath = outputPath;
                obj.postNetwork = fullfile(outputPath,'postsynaptic_network');
                obj.preNetwork = fullfile(outputPath,'presynaptic_network');
                obj.correlationFile = fullfile(obj.preNetwork,'correlations.csv');
                obj.savePath = fullfile(obj.outputPath,'LFPy');
                obj.spikingFile = fullfile(outputPath,'presynaptic_network','spikeTimes.csv');
                mkdir(obj.outputPath);
                mkdir(obj.postNetwork);
                mkdir(obj.preNetwork);
                mkdir(obj.savePath);
            end
        end

        function obj = setsynapsecount(obj,N)
            obj.synapseCount = N;
        end

        function obj = setSpikingFile(obj,filename)
            obj.spikingFile = filename;
        end

        function obj = setSavePath(obj,pathName)
            obj.savePath = pathName;
        end

        function obj = setFiringRate(obj,lambdaE,lambdaI)
            obj.eFiringRate = lambdaE;
            obj.iFiringRate = lambdaI;
        end

        function network = setCorrelationFile(network,file)
            network.correlationFile = file;
        end

        function obj = setPostNetwork(obj,pathName)
            obj.postNetwork = pathName;
        end

        function obj = addActiveChannels(obj,toAdd)
            if(nargin<2)
                obj.activeConductances = true;
            else
                obj.activeConductances = toAdd;
            end
        end

        function obj = changeInhibitorySynapseParams(obj,x)
            % Changes physiology of GABA receptors based on
            % input vector x. x must be of the form
            %    x = [tau_decay,tau_scale_factor,amplitude_scale_factor]
            % (Scaling factors must be less than 10).
            obj.inSynParamChanges = x(:)'*[100;10;1];
        end

        function file = getCorrelationFile(network)
            file = network.correlationFile;
        end

        function N = getsynapsecount(obj)
            N = obj.synapseCount;
        end

        function [time,V,dipoles] = importSimulationResults(network)
            load(fullfile(network.savePath,'simulation_data.mat'));
            time = time(2:end);
            V = V(2:end,:);
            dipoles = dipoles(2:end,:,:);
        end

        function obj = initialize_postsynaptic_network(obj,neuronCount,randMorphos)
        % Initalizes postsynapse network by drawing neuron morphologies and creating empty
        % neuron classes to be filled upon simulation.
           if(isempty(obj.postNetwork))
                error('Directory for postsynaptic network (obj.postNetwork) must be specified.');
            end
            if(exist(obj.postNetwork)==0)
                mkdir(obj.postNetwork);
            end

            % Get distribution of mTypes
            synapseDist = obj.dendriteLengths*1.15;
            mTypes = unique(obj.neuronTypes);
            mTypeCount = accumarray(findgroups(obj.neuronTypes),obj.nrnAbundance);
            mTypeSynapseCount = splitapply(@mean,synapseDist(:),findgroups(obj.neuronTypes));

            % Draw neuron morphologies for network
            if(nargin<3)
                if(isempty(neuronCount))
                    neuronCount = length(mTypes);
                    randMorphos = 1:length(mTypes);
                else
                    randMorphos = pdfSample(mTypeCount,neuronCount);
                end
            end
            obj.morphoCounts = histcounts(randMorphos,0.5:1:length(mTypes)+0.5);
            obj.synapseCount = ceil(sum(obj.morphoCounts(:).*mTypeSynapseCount(:)));
            obj.neurons = struct.empty(neuronCount,0);
            counter = 1;

            % Save morphologies of post synaptic network
            filename3 = fullfile(obj.postNetwork,['mTypes.txt']);
            fid3 = fopen(filename3,'w');
            for i = 1:length(mTypes)
                for k = 1:obj.morphoCounts(i)
                    mPath = fullfile(obj.morphologyPath,[mTypes{i} '.swc']);
                    if(strcmp(filesep,'\'))
                        mPath = strrep(mPath,'\','/');
                    end
                    fprintf(fid3,'%s\n',mPath);
                    simID = ['cell' sprintf('%05d',counter)];
                    obj.neurons(counter).name = simID;
                    obj.neurons(counter).mType = mTypes{i};
                    counter = counter+1;
                end
            end
            fclose(fid3);
        end

        function obj = initialize_presynaptic_network(obj,m,tmax);
        % Generate presynaptic spike times and place synapse locations on postsynaptic network.
        % Postsynapse network must first be initalized.
            obj.branchNo = m;
            if(isempty(obj.synapseCount))
                error('Must initialize postsynaptic network first.')
            end
            if(nargin==3)
                obj.tmax = tmax;
            end
            if(isempty(obj.tmax))
                error('Tmax (ms) must be specified.');
            else
                tmax = obj.tmax;
            end
            if(isempty(obj.preNetwork))
                error('File for presynaptic network (obj.preNetwork) must be specified.');
            end
            if(exist(fileparts(obj.preNetwork))==0)
                mkdir(fileparts(obj.preNetwork));
            else
                filename = obj.spikingFile;
                if(exist(filename))
                    disp('Spike times already exist');
                    return;
                end
            end

            [neuronIDs,spikeTime,ei,B,t] = obj.simulatespikes(m);

            network_simulation_beluga.save_presynaptic_network(neuronIDs,spikeTime,ei,obj.synapseCount,obj.spikingFile);
        end

        function obj = form_connections(obj,coordination_index)
        % Places synapses such that correlated synapses are placed in dendrite segments with similar directions relative to the soma

            if(nargin<2)
                coordination_index = 0;
            end

            % Get IDs for multiple synapses/connection (ensure they end up on the same postsyanptic neurons)
            if(exist(fullfile(obj.preNetwork,'multisynapse_IDs.csv')))
                multID = csvread(fullfile(obj.preNetwork,'multisynapse_IDs.csv'));
            else
                multID = -ones(obj.synapseCount,1);
            end
            mutliK = obj.eiFraction*obj.eMulti+(1-obj.eiFraction)*obj.iMulti;
            multiN = ceil(min(sum(multID==-1),obj.synapseCount)/mutliK);
            parentSynapses = find(multID==-1);

            % Load neuron segment locations projected onto sphere
            if(~exist(obj.mTypeSegmentationData));
                wd = mfilename('fullpath');
                error(['Property *resourceFolder* does not point to data and needs to be changed on line 15 of file ' wd '. If you have not downloaded the data, it is accesible via the link in the README.']);
            end
            load(obj.mTypeSegmentationData);

            % Embed dendrites of each postsyanptic neuron onto sphere
            dendriteEmbedding = cell(length(obj.neurons),1);
            dendriteAreas = cell(length(obj.neurons),1);
            neuronArea = zeros(length(obj.neurons),1);
            for i = 1:length(obj.neurons)
                mData = nrnSegs.(obj.neurons(i).mType);
                x = mean(mData.x,2);
                y = mean(mData.y,2);
                z = mean(mData.z,2);
                D = vecnorm([x,y,z],2,2);
                [theta,phi] = cart2sph(x./D,y./D,z./D);
                dendriteEmbedding{i} = [phi(:),theta(:)];
                dendriteAreaCDF{i} = cumsum(mData.area)/sum(mData.area);
                neuronArea(i) = sum(mData.area);
            end

            % If the syanpses are not placed randomly, perform UMAP on spike train correlation matrix to embed functionally similar presynaptic neurons close to each other
            % Location of synaspes are perturbed away from optimal positions by an amount determined by coordination_index between 0 and 1.
            if(coordination_index>0)
                [synapseEmbedding,idPreSyn] = network_simulation_beluga.loadSynapseLocations(fullfile(obj.preNetwork,'UMAP_embedding.csv'),obj.synapseCount);
                synapseEmbedding = network_simulation_beluga.perturbSynapsePositions(synapseEmbedding,1-coordination_index);
                [theta,phi] = cart2sph(synapseEmbedding(:,1),synapseEmbedding(:,2),synapseEmbedding(:,3));
            else
                M = length(multID);
                theta = asin(2*rand(M,1)-1);
                phi = 2*pi*rand(M,1);
                idPreSyn = 1:M;
            end
            synapseEmbedding = [phi(:),theta(:)];

            % Assign presynaptic neurons to postsynaptic neurons based on total surface area
            % postNeuronList = zeros(obj.synapseCount,1);
            if(length(obj.neurons)>1)
                saCDF = cumsum(neuronArea)/sum(neuronArea);
                postNeuronList = interp1(saCDF,1:length(saCDF),rand(obj.synapseCount,1),'next','extrap');
            else
                postNeuronList = ones(obj.synapseCount,1);
            end

            % Greedy algorithm. Draw dendrites segments and select closest synapse
            segmentList = zeros(obj.synapseCount,1);
            synapseList = zeros(obj.synapseCount,1);
            remainingSyns = true(length(parentSynapses),1);
            multiSynapses = cell(length(obj.neurons),1);
            totalSynapseCount = 0;
            count = 1;
            N = length(parentSynapses);
            while( (count <= N) && (totalSynapseCount <= obj.synapseCount) )
                % Choose random dendrite segment based on surface area
                postNeuron = postNeuronList(count);
                segment = interp1(dendriteAreaCDF{postNeuron},1:length(dendriteAreaCDF{postNeuron}),rand,'next','extrap');

                % Choose presynaptic neuron
                if(any(isnan(dendriteEmbedding{postNeuron}(segment,:))))
                    % If segment is soma, choose a synapse uniformly
                    iClosest = randsample(find(remainingSyns),1);
                else
                    % Find presynaptic neuron closest to segment
                    D = 1-haversine_distance(dendriteEmbedding{postNeuron}(segment,:),synapseEmbedding(parentSynapses,:));
                    [~,iClosest] = max(D.*remainingSyns);
                end
                remainingSyns(iClosest) = false;
                parentSynapseID = idPreSyn(parentSynapses(iClosest));

                % Form synapse
                postNeuronList(count) = postNeuron;
                segmentList(count) = segment;
                synapseList(count) = parentSynapseID;
                count = count + 1;

                % Find all other terminals from presynaptic neuron
                multiSynapses{postNeuron} = [multiSynapses{postNeuron};find(multID(:)==parentSynapseID)];
                totalSynapseCount = length(cat(1,multiSynapses{:}))+count;
            end

            if(totalSynapseCount<obj.synapseCount)
                disp('Warning: fewer presyanptic synapses than expected. Reducing syanpse count. Presyanptic network size too small?');
                postNeuronList = postNeuronList(1:count-1);
                segmentList = segmentList(1:count-1);
                synapseList = synapseList(1:count-1);
                obj.synapseCount = totalSynapseCount;
            end

            % For each postsynaptic neuron, draw segments and add multisynapses
            for postNeuron = 1:length(obj.neurons)
                % Total number of multisynapses
                multiSyns = multiSynapses{postNeuron};
                allSegments = interp1(dendriteAreaCDF{postNeuron},1:length(dendriteAreaCDF{postNeuron}),rand(1,length(multiSyns)),'next','extrap');
                remainingSyns = true(length(multiSyns),1);
                for segment = allSegments
                    % Choose presynaptic neuron
                    if(any(isnan(dendriteEmbedding{postNeuron}(segment,:))))
                        % If segment is soma, choose a synapse uniformly
                        iClosest = randsample(find(remainingSyns),1);
                    else
                        % Find presynaptic neuron closest to segment
                        D = 1-haversine_distance(dendriteEmbedding{postNeuron}(segment,:),synapseEmbedding(multiSyns,:));
                        [~,iClosest] = max(D.*remainingSyns);
                    end
                    remainingSyns(iClosest) = false;
                    parentSynapseID = idPreSyn(multiSyns(iClosest));

                    % Form synapse
                    postNeuronList(count) = postNeuron;
                    segmentList(count) = segment;
                    synapseList(count) = parentSynapseID;
                    count = count + 1;
                end
            end

            % Save connections in seperate csv file
            data = [postNeuronList(:),segmentList(:),synapseList(:)];
            dlmwrite(fullfile(obj.postNetwork,'connections.csv'),data,'delimiter',',','precision','%d');
        end

        function obj = simulate(obj)
            if(isempty(obj.savePath))
                error('Output directory (obj.savePath) must be specified.');
            end
            if(exist(fileparts(obj.savePath))==0)
                mkdir(obj.savePath);
            end
            if(exist(obj.savePath)==0)
                mkdir(obj.savePath)
            end

            if(strcmp(filesep,'\'))
                savePath = strrep(obj.savePath,'\','/');
                functionFolder = strrep(obj.functionFolder,'\','/');
                postNetwork = strrep(obj.postNetwork,'\','/');
            end

            if(~exist(fullfile(obj.morphologyPath,[obj.neuronTypes{1} '.swc'])))
                wd = mfilename('fullpath');
                error(['Property *resourceFolder* does not point to data and needs to be changed on line 15 of file ' wd '. If you have not downloaded the data, it is accesible via the link in the README.']);
            end

            params = strjoin({postNetwork,obj.spikingFile,savePath,int2str(obj.tmax+100),char(string(obj.activeConductances)),num2str(obj.inSynParamChanges)});

            % Convert connections to JSON files
            pySimulate(postNetwork,fullfile(functionFolder,'prep_simulations.py'));

            % Simulate dipoles with LFPy (main simulation)
            pySimulate(params,fullfile(functionFolder,'compute_network_dipoles_beluga.py'));

            % Convert NPY files to mat files
            pySimulate(obj.savePath,fullfile(functionFolder,'npy2mat.py'));

            obj.save();
        end

        function obj = compute_presynaptic_correlations(obj,spikingFile)
            % Get pairwise correlation matrix among all presynaptic
            % neurons based on spike trains (sparse)
            obj.correlationFile = fullfile(obj.preNetwork,'correlations.csv');
            if(exist(obj.correlationFile))
                disp('Correlation file already exists.');
                return;
            end
            if(nargin<2)
                spikingFile = obj.spikingFile;
            end
            [err,prints] = system([obj.correlationFunction ' ' spikingFile ' ' obj.correlationFile]);
            if(err~=0)
                error(sprintf('%d: %s',err,prints));
            end
        end

        function obj = embed_presyanptic_neurons(obj)
            file1 = obj.getCorrelationFile;
            file2 = fullfile(obj.preNetwork,'UMAP_embedding.csv');
            system(['python ' obj.embeddingFunction ' ' file1 ' ' file2]);
        end

        function save(obj,filename)
            if(nargin<2)
                filename = 'data.mat';
            end
            network = obj;
            save(fullfile(obj.outputPath,filename),'network');
        end

        function network2 = copy_network(obj,newPath)
            network3 = network_simulation_beluga(newPath);
            network2 = obj;

            network2.morphologyPath = network3.morphologyPath;
            network2.simulateFunction = network3.simulateFunction;
            network2.correlationFunction = network3.correlationFunction;
            network2.neuronTypes = network3.neuronTypes;
            network2.embeddingFunction = network3.embeddingFunction;

            network2.outputPath = fullfile(newPath);
            network2.preNetwork = fullfile(newPath,'presynaptic_network');
            network2.postNetwork = fullfile(newPath,'postsynaptic_network');
            try
                for i = 1:length(network2.neurons)
                    [~,cellNo] = fileparts(network2.neurons(1).sim);
                    network2.neurons(1).sim = fullfile(newPath,'LFPy',[cellNo '.csv']);
                end
            catch
                disp('No simulation files to copy.');
            end

            copyfile(obj.postNetwork,network2.postNetwork);
            copyfile(obj.spikingFile,fullfile(network2.preNetwork,'spikeTimes.csv'));
            try
                copyfile(fullfile(obj.preNetwork,'UMAP_embedding.csv'),fullfile(network2.preNetwork,'UMAP_embedding.csv'));
                network2.correlationFile = obj.correlationFile;
                fid = fopen(fullfile(network2.preNetwork,'correlations(shortcut).txt'),'w');
                txt = ['Network is using the correlation file at location ' strrep(obj.correlationFile,'\','\\')];
                fprintf(fid,txt);
                fclose(fid);
            catch err
                disp(['Caught exception: ' err.identifier])
            end
            F = dir(obj.savePath);
            F = F(3:end);
            for j = 1:length(F)
                copyfile(fullfile(obj.savePath,F(j).name),fullfile(network2.savePath,F(j).name));
            end
            network2.save();
        end

        function [neuronIDs,spikeTime,ei,B,t] = simulatespikes(obj,m)
            if(m>1)
                error('Branching index must be less than 1.')
            end
            tmax = obj.tmax*1e-3+2.1;
            M = obj.synapseCount;
            dt = 4e-3;
            t = 0:dt:tmax;
            N = length(t);

            ME = floor(obj.eiFraction*M);
            MI = M-ME;

            k = 4;
            h = obj.eFiringRate*dt*ME*(1-m);

            % Generate mean firing rate using critical branching process
            B = zeros(N,1);
            exN = poissrnd(h,N,1);
            for i = 1:N-1
                if(B(i))
                    count = binornd(B(i)*k,m/k);
                else
                    count = 0;
                end
                B(i+1) = count+exN(i);
            end
            B(t<2+dt) = [];
            t(t<2+dt) = []; t = t-2;

            numspikes = sum(B);
            spikeTime = zeros(5*numspikes,1);
            neuronIDs = zeros(5*numspikes,1);
            eiRatio = obj.iFiringRate/obj.eFiringRate;
            j = 0;
            for i = 1:length(t)
                nEx = poissrnd(B(i));
                spikeTime(j+1:j+nEx) = t(i) + rand(1,nEx)*dt - dt/2;
                neuronIDs(j+1:j+nEx) = randperm(ME,nEx);
                j = j+nEx;

                nIn = poissrnd(B(i)/eiRatio);
                spikeTime(j+1:j+nIn) = t(i) + rand(1,nIn)*dt - dt/2;
                neuronIDs(j+1:j+nIn) = ME+randperm(MI,nIn);
                j = j+nIn;
            end
            spikeTime(j+1:end) = [];
            neuronIDs(j+1:end) = [];
            ei = [zeros(1,ME),ones(1,MI)];
            [neuronIDs,I] = sort(neuronIDs);
            spikeTime = 1e3*spikeTime(I);
        end

        function [ids,ts,ei] = simulatespikes_critplane(obj,N,tmax)

            tmax = tmax*1e-3;

            % Number of E and I synapses
            N_ex_syn = floor(obj.eiFraction*N);
            N_in_syn = N-N_ex_syn;

            % Multisynapse count (predicted by Markram et al., Cell 2015)
            N_ex_pre = ceil(N_ex_syn/obj.eMulti);
            N_in_pre = ceil(N_in_syn/obj.iMulti);

            % Initialize
            ei = [zeros(1,N_ex_syn),ones(1,N_in_syn)];
            idsE = -ones(ceil(2*N_ex_pre*tmax),1);
            tsE = -ones(ceil(2*N_ex_pre*tmax),1);

            % Network topology
            elevation_E = rand(N_ex_pre,1);
            azimuth_E = rand(N_ex_pre,1);
            dist_metric = @(x,y)vecnorm(x-y,2,2);
            nNeigh = 10;
            C = zeros(N_ex_pre*nNeigh,3);
            D1 = exprnd(0.02,nNeigh,N_ex_pre);
            for i = 1:N_ex_pre
                D0 = dist_metric([elevation_E(i),azimuth_E(i)],[elevation_E,azimuth_E]);
                D0(1) = Inf;
                idcs = 1:length(D0);
                for j = 1:nNeigh
                    [~,jj] = min(abs(D0(idcs)-D1(j,i)));
                    I(j) = idcs(jj);
                    idcs(jj) = [];
                end
                % [~,I] = sort(D0,'ascend'); I = I(2:nNeigh+1);
                C(nNeigh*(i-1)+1:nNeigh*i,1) = i*ones(nNeigh,1);
                C(nNeigh*(i-1)+1:nNeigh*i,2) = I;
            end
            C(:,3) = obj.branchNo/nNeigh;

            % Look up table for node indices
            for i = 1:N_ex_pre
                CLUT{i} = find(C(:,1)==i);
            end

            % Simulation parameters
            dt = 4e-3;
            t = 0:dt:tmax;
            tN = length(t);

            % Random external input
            nTrans = poissrnd(obj.eFiringRate*dt*N_ex_pre);
            idsE(1:nTrans) = randperm(N_ex_pre,nTrans);
            post = idsE(1:nTrans);
            k = nTrans;
            tsE(1:nTrans) = t(1)*ones(nTrans,1);
            count = nTrans;

            exN = poissrnd(obj.eFiringRate*dt*N_ex_pre*(1-obj.branchNo),tN,1);

            for i = 2:tN-1
                % Get spiking cells at previous time point
                % preI = unique(idsE(count-nTrans+1:count));
                preI = idsE(count-nTrans+1:count);
                jj = cat(1,CLUT{preI});

                % For each spiking cell, find neighbours...
                postI = C(jj,2);
                % ... and start flipping coins
                iTrans = find(rand(length(jj),1)<C(jj,3));
                postI = postI(iTrans);
                nTrans = length(postI);

                % Add propogated spikes
                idsE(count+1:count+nTrans) = postI;
                tsE(count+1:count+nTrans) = t(i)+dt*rand(nTrans,1);;
                count = count+nTrans;

                % Add external noise
                idsE(count+1:count+exN(i)) = randperm(N_ex_pre,exN(i));
                tsE(count+1:count+exN(i)) = t(i)+dt*rand(exN(i),1);
                count = count+exN(i);
                nTrans = nTrans + exN(i);
            end
            tsE(count:end) = [];
            idsE(count:end) = [];

            % Look up table for excitatory spike times
            speLUT = cell(N_ex_pre,1);
            for j = 1:N_ex_pre
                speLUT{j} = find(idsE==j);
            end

            % Inhibitory neurons follow nearby excitatory neurons
            elevation_I = rand(N_in_pre,1);
            azimuth_I = rand(N_in_pre,1);
            idsI = -ones(ceil(10*N_in_pre*tmax),1);
            tsI = -ones(ceil(10*N_in_pre*tmax),1);
            count = 0;

            for i = 1:N_in_pre
                D0 = dist_metric([elevation_I(i),azimuth_I(i)],[elevation_E,azimuth_E]);
                [~,I] = sort(D0,'ascend');
                k = poissrnd(obj.iFiringRate/obj.eFiringRate);
                D1 = exprnd(0.02,k,1);
                D0(1) = Inf;
                idcs0 = 1:length(D0);
                for j = 1:k
                    % Choose excitatory neuron nearby
                    [~,jj] = min(abs(D0(idcs0)-D1(j)));
                    idcs = speLUT{idcs0(jj)};
                    idcs0(jj) = [];

                    nTrans = length(idcs);
                    idsI(count+1:count+nTrans) = i+N_ex_syn;,

                    % Make some spikes random depending on branchNo
                    idcs2 = rand(nTrans,1)<obj.branchNo;
                    t0 = tsE(idcs).*idcs2+tmax*rand(nTrans,1).*(1-idcs2);

                    % Add spikes to inhibitory neuron
                    tsI(count+1:count+nTrans) = t0+dt*rand(nTrans,1);
                    count = count+nTrans;
                end
            end
            tsI(count:end) = [];
            idsI(count:end) = [];

            % Look up table for inhibitory spike times
            spiLUT = cell(N_in_pre,1);
            for j = 1:N_in_pre
                spiLUT{j} = find(idsI==(j+N_ex_syn));
            end

            % Define parent synapses of which multisyanpses are copies
            % If parent is -1, it is the original synapse.
            parents = -ones(N,1);

            % Final excitatory synapses
            M = randi(N_ex_pre,N_ex_syn-N_ex_pre,1);
            for j = 1:length(M)
                % Get random presyanptic neuron
                idcs = speLUT{M(j)};

                % Truncated normal distribution centered on "parent" synapse
                e = truncate(makedist('Normal','mu',elevation_E(M(j)),'sigma',0.02),0,1);
                a = truncate(makedist('Normal','mu',azimuth_E(M(j)),'sigma',0.02),0,1);

                elevation_E(j+N_ex_pre) = e.random; % Set synapse location nearby
                azimuth_E(j+N_ex_pre) = a.random; % Set synapse location nearby
                parents(j+N_ex_pre) = M(j);

                tsE = [tsE;tsE(idcs)]; % Copy spike times
                idsE = [idsE;(j+N_ex_pre)*ones(length(idcs),1)]; % New synapse ID
            end

            % Final inhibitory synapses
            M = randi(N_in_pre,N_in_syn-N_in_pre,1);
            for j = 1:length(M)
                % Get random presyanptic neuron
                idcs = spiLUT{M(j)};

                % Truncated normal distribution centered on "parent" synapse
                e = truncate(makedist('Normal','mu',elevation_I(M(j)),'sigma',0.05),0,1);
                a = truncate(makedist('Normal','mu',azimuth_I(M(j)),'sigma',0.05),0,1);

                elevation_I(j+N_in_pre) = e.random; % Set synapse location nearby
                azimuth_I(j+N_in_pre) = a.random; % Set synapse location nearby
                parents(j+N_in_pre+N_ex_syn) = M(j)+N_ex_syn;

                tsI = [tsI;tsI(idcs)]; % Copy spike times
                idsI = [idsI;(j+N_in_pre+N_ex_syn)*ones(length(idcs),1)]; % New synapse ID
            end

            % Combine
            elevation = [elevation_E;elevation_I];
            azimuth = [azimuth_E;azimuth_I];
            ids = [idsE;idsI];
            ts = [tsE;tsI]*1e3;

            csvwrite(fullfile(obj.preNetwork,'connections.csv'),C);
            csvwrite(fullfile(obj.preNetwork,'locations.csv'),[elevation(:),azimuth(:)]);
            csvwrite(fullfile(obj.preNetwork,'multisynapse_IDs.csv'),parents(:));
            network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,obj.spikingFile)
        end

    end

    methods (Static)

        function [embedding,idPreSyn] = loadSynapseLocations(umapFile,N)

            data = csvread(umapFile);
            [x,y,z] = sph2cart(data(:,2),data(:,3),data(:,3)*0+1);
            idPreSyn = data(:,1)+1;
            embedding = [x(:),y(:),z(:)];
            M = size(embedding,1);

            if(M<N)
                % Place all the left over synapses, i.e. those that were not correlated, randomly on the sphere
                K = N-M;
                elevation = asin(2*rand(K,1)-1);
                azimuth = 2*pi*rand(K,1);
                [x,y,z] = sph2cart(azimuth,elevation,1);
                uncorrelated = [x,y,z];
                embedding = [embedding;uncorrelated];
                idPreSyn0 = setdiff(1:N,idPreSyn);
                idPreSyn = [idPreSyn;idPreSyn0(:)];
            end
        end

        function save_presynaptic_network(neuronIDs,spikeTime,ei,synapseCount,filename)
            % Save spiketimes generated by simulatespikes as a csv file for later use
            fid = fopen(filename,'w');
            [uM,~,nTemp] = unique(neuronIDs(:));
            I = accumarray(nTemp,1:numel(nTemp),[],@(x){sort(x)});
            UC = 0;
            for i = 1:synapseCount
                if(i-UC>length(uM))
                    ts = [];
                    ei(i) = rand>0.85;
                else
                    if(uM(i-UC)==i)
                        ts = spikeTime(I{i-UC});
                    else
                        ts = [];
                        UC = UC+1;
                    end
                end
                if(ei(i)==0)
                    eiType = 'e';
                else
                    eiType = 'i';
                end
                s = sprintf('%0.1f,',sort(ts));
                fprintf(fid,'%s,%s\n',eiType,s(1:end-1));
                % fprintf(fid,'%s\n',s(1:end));
            end
            fclose(fid);
        end

         function [f,P] = expectedEEGspectrum(dipoles,sa,idcs)
            for j = 1:size(dipoles,3)
                temp = resample(dipoles(2:end,:,j),2e3,16e3);
                if(j==1)
                    dp = zeros(size(temp,1),3,size(dipoles,3));
                end
                dp(:,:,j) = temp;
            end
            q = gpuArray(dp);

            for j = 1:length(idcs)
                eeg = detrend(network_simulation_beluga.getEEG(q,sa,idcs(j)),'constant');
                temp = mypmtm(eeg,2e3,2);
               if(j==1)
                    P = zeros(size(temp));
                end
                P = P+temp;
            end
            P = P/length(idcs)*pi/2;
            f = 0.5:0.5:1e3;
        end

        function W = getEEG(Q,sa,idx)
            % Lead field for Cz recording
            L0 = squeeze(sa.cortex75K.V_fem(49,idx,:))'; % czIDX = 49

            % Orient dipole "up" direction to be normal to cortical surface
            vz = sa.cortex75K.normals(idx,:);
            [vx,vy] = getOrthBasis(vz);
            q = (Q(:,1,:).*vx+Q(:,2,:).*vy+Q(:,3,:).*vz)*1e-12; % nAum -> mAm

            W = squeeze(sum(L0.*q,2)*1e6); % V -> uV
        end

        function [sa,X] = getHeadModel()
            try
                load(fullfile(network_simulation_beluga.resourceFolder,'anatomy_nyhead_model.mat'));
            catch
                wd = mfilename('fullpath');
                error(['Could not read file. Property *resourceFolder* does not point to data and needs to be changed on line 15 of file ' wd '. If you have not downloaded the data, it is accesible via the link in the README.']);
            end
            X = struct();
            X.vertices = sa.cortex75K.vc;
            X.faces= sa.cortex75K.tri;
        end

        function [neuronIDs,spikeTime,ei] = getprenetwork(spikingFile)
            T = dlmread(spikingFile,',',0,1);
            neuronIDs = repmat((1:size(T,1)),[size(T,2),1])';
            spikeTime = T(T~=0);
            neuronIDs = neuronIDs(T~=0);

            fid = fopen(spikingFile);
            T = textscan(fid,'%s');
            ei = cellfun(@(x)x(1),T{1})=='i';
            fclose(fid);
        end

        function Xnew = perturbSynapsePositions(X,S,isCart)

            if(nargin<3 || isCart)
                isCart = true;
                [theta0,phi0] = cart2sph(X(:,1),X(:,2),X(:,3));
                phi0 = phi0+pi/2;
            elseif(~isCart)
                theta0 = X(:,1);
                phi0 = X(:,2)+pi/2;
            end

            % Set a random distance and random bearing (distance propotional to S)
            sz = size(phi0);
            B = pi*rand(sz);
            u = 1-2*S*rand(sz);
            idcs = (rand(sz)<0.5);
            d = acos(u).*(1-2*idcs)+2*pi*idcs;

            % Compute final locations of synapses
            % d = linspace(0,2*pi,100);
            phi2 = acos(cos(d).*cos(phi0)+sin(d).*sin(phi0).*cos(B));
            A = acos((cos(d)-cos(phi2).*cos(phi0))./(sin(phi2).*sin(phi0)));
            A(d>pi) = 2*pi-A(d>pi);
            theta2 = theta0+real(A);

            if(isCart)
                [x,y,z] = sph2cart(theta2,phi2-pi/2,1);
                Xnew = [x,y,z];
            else
                Xnew = [theta2,phi2-pi/2];
            end
        end

        function d = haversine_distance2(y,x)
            d = hav(y(:,1)-x(:,1))+(1-hav(x(:,1)-y(:,1))-hav(y(:,1)+x(:,1))).*hav(y(:,2)-x(:,2));
        end
    end
end
function prints = pySimulate(params,pyFun)
    [err,prints] = system(['python "' pyFun '" ' params]);
    if(err)
        error(prints);
    end
end
function samples = pdfSample(empiricalDist,n)
    occuranceCDF = cumsum(empiricalDist)/sum(empiricalDist);
    samples = interp1(occuranceCDF,1:length(occuranceCDF),rand(n,1),'next','extrap');
end
function d = haversine_distance(p1,x)
    d = hav(p1(1)-x(:,1))+(1-hav(x(:,1)-p1(1))-hav(p1(1)+x(:,1))).*hav(p1(2)-x(:,2));
end
function d = hav(x)
    d = (1-cos(x))/2;
end
