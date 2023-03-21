classdef simneuron

    properties
        sim
        morph
        syn
        savePath
        name
    end

    properties (Access = private)
        morphDataPath = 'E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\cortical_column_Hagen\segment_areas.mat'
    end

    methods
        function obj = simneuron(name,savePath,neuronID);
            obj.sim = struct('V_soma',[],'time',[],'dipole',[],'fs',[]);
            obj.morph = struct('x',[],'neuronID','');
            obj.syn = struct('x',[],'x_norm',[],'n',[],'IDs',[],'preSyns','');
            if(nargin>=1)
                obj.name = name;
            end
            if(nargin>=2)
                obj = obj.load(savePath);
            end
            if(nargin>=3)
                obj = obj.addMorphology(neuronID);
            end
        end
        function obj = load(obj,savePath)

            data = load(savePath);
            obj.savePath = savePath;

            obj.sim = struct('V_soma',[],'time',[],'dipole',[],'fs',[]);
            obj.sim.V_soma = data.V_soma(:);
            obj.sim.time = data.t(:);
            obj.sim.fs = 1e3./mean(diff(data.t));
            obj.sim.dipole = [data.Q(1,:)',-data.Q(3,:)',data.Q(2,:)'];

            synapseFile = data.synapseFile;
            synInfo = csvread(synapseFile);

            obj.syn = struct('x',[],'x_norm',[],'n',size(synInfo,1),'IDs',synInfo(:,1),'preSyns',synInfo(:,2));

            obj.morph = struct('x',[],'neuronID','','area','');
        end
        function obj = addMorphology(obj,neuronID)
            load(obj.morphDataPath)
            obj.morph.x = nrnSegSA.(neuronID).pos;
            obj.morph.neuronID = neuronID;
            obj.morph.area = nrnSegSA.(neuronID).area;
        end
        function obj = addSynapses(obj,synapseFile)
            synInfo = csvread(synapseFile);
            obj.syn = struct('x',[],'x_norm',[],'n',size(synInfo,1),'IDs',synInfo(:,1),'preSyns',synInfo(:,2));
            if(isempty(obj.morph.x))
                error('Mophology must be initialized prior to synapse addition.');
            else
                obj = obj.addSynapseLocations;
            end
        end
        function obj = addSynapseLocations(obj)
            obj.syn.x = obj.morph.x(obj.syn.IDs,:);
            obj.syn.x_norm = vec2norm(obj.syn.x);
        end
        function similarityIndex = synDist(obj,neuron2)
            iU = intersect(obj.syn.preSyns,neuron2.syn.preSyns);
            iT = union(obj.syn.preSyns,neuron2.syn.preSyns);
            d = zeros(length(iU),1);
            for i = 1:length(iU)
                idx1 = find(obj.syn.preSyns==iU(i));
                idx2 = find(neuron2.syn.preSyns==iU(i));
                d(i) = dot(obj.syn.x_norm(idx1,:),neuron2.syn.x_norm(idx2,:));
            end

            d(isnan(d))=1;
            similarityIndex = sum(d)/length(iT);
        end
        function uP = getUnitarySpectrum(obj,sa)
            idcs3 = sa.cortex2K.in_from_cortex75K;
            for i = 1:length(idcs3)
                eeg = obj.getEEG(sa,idcs3(i));
            end
        end
        function eeg = getEEG(obj,sa,idx)
            eeg = getEEG(obj.syn.dipole,sa,idx);
        end

    end
end

function X_norm = vec2norm(X)
    X_norm = X./vecnorm(X,2,2);
end
