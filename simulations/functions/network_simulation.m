classdef network_simulation

    properties
        tmax
        branchingIdx
        neurons
        outputPath
        postNetwork
        preNetwork
        spikingFile
        savePath
        results = struct('t',[],'dipoles',[],'Vmem',[]);
    end

    properties (Constant)
        resourceFolder = 'E:/Research_Projects/004_Propofol/data/resources';
        functionFolder = fileparts(mfilename('fullpath'));
    end

    properties (Access = private)
        morphologyPath = fullfile(network_simulation.resourceFolder,'cortical_column_Hagen','swc');
        mTypeSegmentationData = fullfile(network_simulation.resourceFolder,'cortical_column_Hagen','segment_areas.mat');
        simulateFunction = fullfile(network_simulation.functionFolder,'compute_network_dipoles.py');
        correlationFunction = fullfile(network_simulation.functionFolder,'compute_correlation_matrix_parallel.exe');
        embeddingFunction = fullfile(network_simulation.functionFolder,'embed_data.py');

        neuronTypes = {'L23E_oi24rpy1';'L23I_oi38lbc1';'L23I_oi38lbc1';'L4E_53rpy1';'L4E_j7_L4stellate';'L4E_j7_L4stellate';'L4I_oi26rbc1';'L4I_oi26rbc1';'L5E_oi15rpy4';'L5E_j4a';'L5I_oi15rbc1';'L5I_oi15rbc1';'L6E_51_2a_CNG';'L6E_oi15rpy4';'L6I_oi15rbc1';'L6I_oi15rbc1'};
        nrnAbundance = [26.8,3.2,4.3,9.5,9.5,9.5,5.6,1.5,4.9,1.3,0.6,0.8,14,4.6,1.9,1.9]/100
        dendriteLengths = [9319,4750,4750,7169,4320,4320,2359,2359,14247,5720,5134,5134,5666,6183,5134,5134];
        samplingFrequency = 16000
        activeConductances = false
        propofol = 0
        synapseCount
        morphoCounts
        correlationFile
        eiFraction = 0.85;
        eFiringRate = 1; % Hz
        iFiringRate = 6.67; % Hz
    end

    methods

        function obj = network_simulation(outputPath)
            warning('off','MATLAB:MKDIR:DirectoryExists')
            if(nargin>0)
                obj.outputPath = outputPath;
                obj.postNetwork = fullfile(outputPath,'postsynaptic_network');
                obj.preNetwork = fullfile(outputPath,'presynaptic_network');
                obj.correlationFile = fullfile(obj.preNetwork,'correlations.csv');
                mkdir(obj.outputPath);
                mkdir(obj.postNetwork);
                mkdir(obj.preNetwork);
                mkdir(fullfile(obj.postNetwork,'connections'));
                mkdir(fullfile(obj.outputPath,'LFPy'));
            end
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
            if(exist(fullfile(obj.postNetwork,'connections'))==0)
                mkdir(fullfile(obj.postNetwork,'connections'));
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
                    fprintf(fid3,'%s\n',fullfile(obj.morphologyPath,[mTypes{i} '.swc']));
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
            obj.branchingIdx = m;
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
                filename = fullfile(obj.preNetwork,'spikeTimes.csv');
                if(exist(filename))
                    disp('Spike times already exist');
                    return;
                end
            end

            [neuronIDs,spikeTime,ei,B,t] = obj.simulatespikes(m);

            network_simulation.save_presynaptic_network(neuronIDs,spikeTime,ei,obj.synapseCount,fullfile(obj.preNetwork,'spikeTimes.csv'));
        end

        function form_connections(obj,coordination_index)
        % Places synapses such that correlated synapses are placed in dendrite segments with similar directions relative to the soma

            if(nargin<2)
                coordination_index = 0;
            end

            % Load neuron segment locations projected onto sphere
            load(obj.mTypeSegmentationData)
            mTypes = unique(obj.neuronTypes);
            mPresent = find(obj.morphoCounts>0);
            mTypes = mTypes(mPresent);
            morphCounts = obj.morphoCounts(mPresent);
            morphCounts = morphCounts(:);

            X = []; Y = []; Z = []; nrnID = []; segID = []; SA = {}; contIdx = {};
            counter = 0;
            counter2 = 1;
            for i = 1:length(mTypes)
                mData = nrnSegs.(mTypes{i});
                N = length(mData.area);
                pos = [mean(mData.x,2),mean(mData.y,2),mean(mData.z,2)];
                segDir = pos./vecnorm(pos,2,2);
                nullIdcs = find(sum(isnan(segDir),2));
                segDir(nullIdcs,:) = segDir(randi(N,length(nullIdcs)),:);
                X = [X;segDir(:,1)];
                Y = [Y;segDir(:,2)];
                Z = [Z;segDir(:,3)];
                segID = [segID;(1:N)'];
                SA{i} = mData.area(:);
                contIdx{i} = (1:N)'+counter;
                counter = counter+N;
                for j = 1:morphCounts(i)
                    nrnID = [nrnID; counter2*ones(N,1)];
                    counter2 = counter2+1;
                end
            end
            dendriteEmbedding = [X,Y,Z];

            % Randomly draw neuron segments from postsynaptic population,
            % proportionally to segment surface area
            mcs = mat2cell(morphCounts,ones(length(contIdx),1))';
            sa = cellfun(@(x,y) repmat(x,[y,1]),SA,mcs,'UniformOutput',false);
            sa = cat(1,sa{:});
            contIdx = cellfun(@(x,y) repmat(x,[y,1]),contIdx,mcs,'UniformOutput',false);
            contIdx = cat(1,contIdx{:});


            saCDF = cumsum(sa)/sum(sa);
            iSegsPooled = interp1(saCDF,1:length(saCDF),rand(obj.synapseCount,1),'next','extrap');
            iSegs = contIdx(iSegsPooled);

            % If the syanpses are not placed randomly, perform UMAP on spike train correlation matrix to embed functionally similar presynaptic neurons close to each other
            % Location of synaspes are perturbed away from optimal positions by an amount determined by coordination_index between 0 and 1.
            if(coordination_index>0)
                synapseEmbedding = obj.embedPresynapticNetwork();
                synapseEmbedding = network_simulation.perturbSynapsePositions(synapseEmbedding,1-coordination_index);

                % Convert to spherical coordinates to compute Haversine distances
                [theta,phi] = cart2sph(synapseEmbedding(:,1),synapseEmbedding(:,2),synapseEmbedding(:,3));
                synapseEmbedding = [phi(:),theta(:)];
                [theta,phi] = cart2sph(dendriteEmbedding(:,1),dendriteEmbedding(:,2),dendriteEmbedding(:,3));
                dendriteEmbedding = [phi(:),theta(:)];

                % Greedy algorithm. For each neuron segment select the presynaptic neuron closest on the sphere embedding
                syn = zeros(size(iSegs));
                remainingSyns = true(size(synapseEmbedding,1),1);
                h = waitbar(0,'Wiring network together...');
                for i = 1:obj.synapseCount
                    waitbar(i/obj.synapseCount,h)
                    D = 1-haversine_distance(dendriteEmbedding(iSegs(i),:),synapseEmbedding);
                    [~,syn(i)] = max(D.*remainingSyns);
                    remainingSyns(syn(i)) = false;
                end
                delete(h)
            else
                syn = randperm(obj.synapseCount);
            end

            % Save connections in seperate csv files for each neuron
            % First column inidicates the segment number and the
            % second column indicates the pre synaptic neuron ID
            data = zeros(obj.synapseCount,3);
            for i = 1:obj.synapseCount
                j = nrnID(iSegsPooled(i));
                seg = segID(iSegs(i));
                pre = syn(i);
                data(i,:) = [j,seg,pre];
            end
            for j = 1:length(obj.neurons)
                simID = obj.neurons(j).name;
                filename = fullfile(obj.postNetwork,'connections',[simID '.csv']);
                fid = fopen(filename,'w');
                idcs = find(data(:,1)==j);
                for i = 1:length(idcs)
                    fprintf(fid,'%s,%s\n',int2str(data(idcs(i),2)),int2str(data(idcs(i),3)));
                end
                fclose(fid);
            end
        end

        function embedding = embedPresynapticNetwork(obj)
            N = obj.synapseCount;
            embedding = zeros(N,3);

            correlationFile = obj.correlationFile;
            if(~exist(correlationFile))
                error('Correlation file must be pre computed using function "compute_presynaptic_correlations"');
            end

            % Run UMAP using python module
            % Haversine metric allows embedding data onto a sphere
            pyFun = obj.embeddingFunction;
            umapFile = fullfile(fileparts(network.getCorrelationFile),'UMAP_embedding.csv');
            if(~exist(umapFile))
                [err,prints] = system(['python "' pyFun '" ' correlationFile ' ' umapFile]);
                if(err)
                    error(prints);
                end
            end
            data = csvread(umapFile);
            [x,y,z] = sph2cart(data(:,2),data(:,3)+pi/2,data(:,3)*0+1);
            embedding = [x(:),y(:),z(:)];
            M = size(embedding,1);

            if(M>N)
                embedding = embedding(1:N,:);
            elseif(M<N)
                % Place all the left over synapses, i.e. those that were not correlated, randomly on the sphere
                K = N-M;
                elevation = asin(2*rand(K,1)-1);
                azimuth = 2*pi*rand(K,1);
                [x,y,z] = sph2cart(azimuth,elevation,1);
                uncorrelated = [x,y,z];
                embedding = [embedding;uncorrelated];
            end
        end

        function obj = compute_presynaptic_correlations(obj,spikingFile)
            % Get pairwise correlation matrix among all presynaptic
            % neurons based on spike trains (sparse)
            if(exist(obj.correlationFile))
                disp('Correlation file already exists.');
                return;
            end
            if(nargin<2)
                spikingFile = fullfile(obj.preNetwork,'spikeTimes.csv');
            end
            [err,prints] = system([obj.correlationFunction ' ' spikingFile ' ' obj.correlationFile]);
            if(err~=0)
                error(sprintf('%d: %s',err,prints));
            end
        end

        function obj = addActiveChannels(obj,toAdd)
            if(nargin<2)
                obj.activeConductances = true;
            else
                obj.activeConductances = toAdd;
            end
        end

        function obj = addPropofol(obj,x)
            if(nargin<2)
                obj.propofol = 1;
            else
                obj.propofol = x;
            end
            % obj.eFiringRate = obj.eFiringRate*0.5;
        end

        function params = getSimulationParams(obj)
            params.postNetwork = obj.postNetwork;
            params.spikingFile = fullfile(obj.preNetwork,'spikeTimes.csv');
            params.outputPath = obj.outputPath;
            params.tmax = obj.tmax+100;
            params.soma = char(string(obj.activeConductances))
            params.propofol = obj.propofol;
            params.simulateFunction = obj.simulateFunction;
        end

        function obj = simulate(obj)
            if(isempty(obj.outputPath))
                error('Output directory (obj.outputPath) must be specified.');
            end
            if(exist(fileparts(obj.outputPath))==0)
                mkdir(obj.outputPath);
            end
            if(exist(fullfile(obj.outputPath,'LFPy'))==0)
                mkdir(fullfile(obj.outputPath,'LFPy'));
            end

            spikingFile = fullfile(obj.preNetwork,'spikeTimes.csv');
            params = strjoin({obj.postNetwork,spikingFile,obj.outputPath,int2str(obj.tmax+100),char(string(obj.activeConductances)),num2str(obj.propofol)});
            prints = pySimulate(params,obj.simulateFunction);
            outputFiles = split(prints,char(10));
            for j = 1:length(obj.neurons)
                id = obj.neurons(j).name;
                idx = find(cellfun(@(x)~isempty(strfind(x,id)),outputFiles));
                obj.neurons(j).sim = outputFiles{idx};
            end
            obj.importResults;
            obj.save()
        end

        function save(obj)
            network = obj;
            save(fullfile(obj.outputPath,'data.mat'),'network');
        end

        function N = getsynapsecount(obj)
            N = obj.synapseCount;
        end

        function obj = setsynapsecount(obj,N)
            obj.synapseCount = N;
        end

        function obj = importResults(obj)
            k = 1;
            n = length(obj.neurons);
            for i = 1:n
                data = csvread(obj.neurons(i).sim,1,0);
                t = data(:,1);
                if(i==1)
                    Q = zeros([size(data,1),3,n]);
                    V = zeros([size(data,1),n]);
                end
                Q(:,:,k) = data(:,2:4);
                V(:,k) = data(:,5);
                k = k+1;
            end
            Q = [Q(:,1,:),-Q(:,3,:),Q(:,2,:)];
            Q = Q(t>=100,:,:);
            V = V(t>=100,:,:);
            t = t(t>=100); t = t-100;
            obj.results.t = t;
            obj.results.dipoles = Q;
            obj.results.Vmem = V;
        end

        function [f,P] = expectedEEGspectrum(obj,neuronIDs,sa,idcs)
            if(isempty(neuronIDs))
                neuronIDs = 1:length(obj.neurons);
            end
            if(nargin<4)
                idcs = sa.cortex2K.in_from_cortex75K;
            end

            for j = 1:length(neuronIDs);
                temp = resample(obj.results.dipoles(2:end,:,neuronIDs(j)),2e3,16e3);
                if(j==1)
                    dp = zeros(size(temp,1),3,length(neuronIDs));
                end
                dp(:,:,j) = temp;
            end
            q = gpuArray(dp);
            t = obj.results.t;

            h = waitbar(0);
            for j = 1:length(idcs)
                waitbar(j/length(idcs),h);
                eeg = network_simulation.getEEG(q,sa,idcs(j));
                temp = mypmtm(eeg,2e3,2);
               if(j==1)
                    P = zeros(size(temp));
                end
                P = P+temp;
            end
            close(h);
            P = P/length(idcs)*pi/2;
            f = 0.5:0.5:1e3;
        end

        function [f,P] = getUnitarySpectrum(obj,neuronIDs,sa,idcs,toSum)
            if(isempty(neuronIDs))
                neuronIDs = 1:length(obj.neurons);
            end
            if(nargin<3)
                idcs = sa.cortex2K.in_from_cortex75K;
            end
            if(nargin==5 && toSum)
                q = sum(obj.results.dipoles(:,:,neuronIDs),3);
            else
                q = obj.results.dipoles(:,:,neuronIDs);
            end
            h = waitbar(0);
            for i = 1:length(idcs)
                waitbar(i/length(idcs),h);
                eeg = network_simulation.getEEG(q,sa,idcs(i));
                [temp,f] = pmtm(eeg,2,[],1e3*16);
                if(i==1);
                    P = zeros(size(temp));
                end
                P = P + temp;
            end
            close(h);
            P = P/length(idcs);
        end

        function [neuronIDs,spikeTime,ei] = getprenetwork(obj,file)
            if(nargin<2)
                file = fullfile(obj.preNetwork,'spikeTimes.csv');
            end
            spikingFile = fullfile(file);
            T = dlmread(spikingFile,',',0,1);
            neuronIDs = repmat((1:size(T,1)),[size(T,2),1])';
            spikeTime = T(T~=0);
            neuronIDs = neuronIDs(T~=0);

            fid = fopen(spikingFile);
            T = textscan(fid,'%s');
            ei = cellfun(@(x)x(1),T{1})=='i';
            fclose(fid);
        end

        function plotSynapses(obj,nrnID,clrs)
            load(obj.mTypeSegmentationData)
            nrn = obj.neurons(nrnID);
            mData = nrnSegs.(nrn.mType);
            cons = csvread(fullfile(obj.postNetwork,'connections',[nrn.name '.csv']));

            pos = [mean(mData.x,2),mean(mData.y,2),mean(mData.z,2)];

            figureNB;
            scatter3(pos(cons(:,1),1),pos(cons(:,1),2),pos(cons(:,1),3),5,clrs(cons(:,2),:),'filled');
            hold on;
            line(mData.x',mData.y',mData.z','color','k')
            set(gca,'DataAspectRatio',[1,1,1]);
            gcaformat
        end

        function obj = setFiringRate(obj,lambdaE,lambdaI)
            obj.eFiringRate = lambdaE;
            obj.iFiringRate = lambdaI;
        end
        function network = setCorrelationFile(network,file)
            network.correlationFile = file;
        end

        function file = getCorrelationFile(network)
            file = network.correlationFile;
        end

        function network2 = copy_network(obj,newPath)
            network3 = network_simulation(newPath);
            network2 = obj;

            network2.morphologyPath = network3.morphologyPath;
            network2.simulateFunction = network3.simulateFunction;
            network2.correlationFunction = network3.correlationFunction;
            network2.mTypeSegmentationData = network3.mTypeSegmentationData;
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
            copyfile(fullfile(obj.preNetwork,'spikeTimes.csv'),fullfile(network2.preNetwork,'spikeTimes.csv'));
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
            F = dir(fullfile(obj.outputPath,'LFPy'));
            F = F(3:end);
            for j = 1:length(F)
                copyfile(fullfile(obj.outputPath,'LFPy',F(j).name),fullfile(network2.outputPath,'LFPy',F(j).name));
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
            B(t<2) = [];
            t(t<2) = []; t = t-2;

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

    end

    methods (Static)

        % function obj = simulate(params)
        %     strParams = strjoin({params.postNetwork,params.spikingFile,params.outputPath,int2str(params.tmax),params.soma,num2str(params.propofol)});
        %     save(fullfile(params.outputPath,'params.mat'));
        %     pySimulate(strParams,params.simulateFunction);
        % end

        % Generate network object from directory
        % function obj = reap_network(folder)
        %     network = network_simulation(folder);
        %     for j = 1:length(obj.neurons)
        %         id = obj.neurons(j).name;
        %         idx = find(cellfun(@(x)~isempty(strfind(x,id)),outputFiles));
        %         obj.neurons(j).sim = outputFiles{idx};
        %     end
        %     obj.importResults;
        %     obj.save()
        % end

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

        function W = getEEG(Q,sa,idx)
            % Lead field for Cz recording
            L0 = squeeze(sa.cortex75K.V_fem(49,idx,:))'; % czIDX = 49

            % Orient dipole "up" direction to be normal to cortical surface
            vz = sa.cortex75K.normals(idx,:);
            [vx,vy] = getOrthBasis(vz);
            q = (Q(:,1,:).*vx+Q(:,2,:).*vy+Q(:,3,:).*vz)*1e-12; % nAum -> mAm

            W = squeeze(sum(L0.*q,2)*1e6); % V -> uV
        end

        function plotspectrum(f,P,varargin)
            idcs = find(and(f>0,f<150));
            plotwitherror(f(idcs),P(idcs,:),'Q',varargin{:});
            set(gca,'xscale','log');
            set(gca,'yscale','log');
            % i = get(gca,'YTick');
            % yls = split(sprintf('10^{%d},',i),',');
            % yticklabels(yls(1:end-1));
            xlim([0.5,150]);
            xticks([1,10,100]);
            xticklabels([1,10,100]);
            xlabel('Frequency (Hz)');
            ylabel(['PSD (' char(956) 'V^2/Hz)'])
            gcaformat;
        end

        function [sa,X] = getHeadModel()
            load(fullfile(network_simulation.resourceFolder,'head_models','sa_nyhead.mat'));
            X = struct();
            X.vertices = sa.cortex75K.vc;
            X.faces= sa.cortex75K.tri;
        end

        function [t,V,Q] = importSimulationFile(filename)
            data = csvread(filename,1,0);
            t = data(:,1);
            Q = data(:,2:4);
            V = data(:,5);
        end

        function Xnew = perturbSynapsePositions(X,S)

            [theta0,phi0] = cart2sph(X(:,1),X(:,2),X(:,3));
            phi0 = phi0+pi/2;

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

            [x,y,z] = sph2cart(theta2,phi2-pi/2,1);
            Xnew = [x,y,z];
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
