EIratio = 10;
L =  10.^linspace(-2,1,7);
for lambda0 = L
    folder = fullfile('E:\Research_Projects\004_Propofol\data\simulations\raw\model_comparison',['EI' int2str(EIratio)],sprintf('frequency_%.2fHz',lambda0));
    [time,dipoles_simulation] = simulate_neuron(folder,lambda0,EIratio);
    [T,dp] = analytical_solution(folder);
    dipoles_computed = interp1(T,dp,time);
    save(fullfile(folder,'model_comparison.mat'),'dipoles_computed','dipoles_simulation','time');
end

function [time,dipoles] = simulate_neuron(folder,lambda0,EIratio)
    % Initialize network
    network = network_simulation_beluga(folder);

    % Initialize post network
    nPostNeurons = 10;
    network = network.initialize_postsynaptic_network(nPostNeurons);

    % Presyanptic network parameters
    nPreNeurons = network.getsynapsecount;
    network.tmax = 2e3; % 2 seconds
    network.branchNo = 0;


    tmax = network.tmax+100;
    nE = floor(network.eiFraction*nPreNeurons);
    nI = nPreNeurons-nE;
    nEx = poissrnd(lambda0*nE*tmax/1e3);
    nIn = poissrnd(EIratio*lambda0*nI*tmax/1e3);
    ei = [zeros(nE,1);ones(nI,1)];
    ids = [randi(nE,nEx,1);nEx+randi(nI,nIn,1)];
    ts = tmax*[rand(nEx,1);rand(nIn,1)];
    network_simulation_beluga.save_presynaptic_network(ids,ts,ei,nPreNeurons,network.spikingFile);

    % Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
    network = network.form_connections(0);

    % Simulate dipoles
    network = network.simulate();

    load(fullfile(folder,'LFPy','simulation_data.mat'));
    time = time(2:end);
    dipoles = dipoles(2:end,:,:);
end
function [T,dp] = analytical_solution(folder)
    load(fullfile(folder,'data.mat'))
    load('E:/Research_Projects/004_Propofol/manuscript/Version3/Data/cortical_column_Hagen/morphology_segmentations.mat','nrnSegs');
    synCons = csvread(fullfile(folder,'postsynaptic_network\connections.csv'));
    [ids,ts,ei] = network_simulation_beluga.getprenetwork(network.spikingFile);

    fid = fopen(fullfile(folder,'postsynaptic_network','mTypes.txt'));
    mTypes = textscan(fid,'%s');
    [~,mTypes] = cellfun(@(x)fileparts(x),mTypes{1},'UniformOutput',false);
    fclose(fid);
    for i = 1:length(mTypes)
        mData = nrnSegs.(mTypes{i});
        pos{i} = [mean(mData.x,2),mean(mData.y,2),mean(mData.z,2)];
    end

    synapses = cell(size(synCons,1),3);
    cellIDs = synCons(:,1);
    for i = 1:size(synCons,1)
        synapses{i,1} = pos{cellIDs(i)}(synCons(i,2),:);
        idcs = find(ids==synCons(i,3));
        synapses{i,2} = ts(idcs);
        synapses{i,3} = ei(synCons(i,3));
    end

    dt = 0.1;
    T = 0:dt:2e3;
    tcov = 0:dt:200;
    covFun = @(x,tau1,tau2,A) A*(exp(-x/tau1)-exp(-x/tau2))/((1-tau2/tau1)*(tau2/tau1).^(-tau2/(tau2-tau1))).*(x>=0);
    icov = covFun(tcov(:)-0.06,28.5,2.4,0.54);
    ecov = covFun(tcov(:)-0.06,1.67,1,1.09);

    u = unique(cellIDs);
    dp = zeros(length(T),3,length(u));

    waitbar(0);
    for i = 1:size(synapses,1)
        waitbar(i/size(synapses,1));
        p = synapses{i,1};
        r = norm(p);
        if(r==0)
            continue;
        end
        p = p/r;
        eventTimes = interp1(T,1:length(T),synapses{i,2},'nearest');
        eventTimes(isnan(eventTimes)) = [];
        events = zeros(size(T));
        events(eventTimes) = 1;
        if(synapses{i,3}==0)
            PSPs = filter(ecov,1,events);
        else
            PSPs = filter(icov,1,events);
        end
        dp(:,:,cellIDs(i)) = dp(:,:,cellIDs(i))+PSPs(:).*p;
    end
    dp = -dp;
end