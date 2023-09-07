function [ids,ts,ei] = resimulate_critplane(network,tmax)


    C = csvread(fullfile(obj.preNetwork,'connections.csv'));
    locations = csvread(fullfile(obj.preNetwork,'locations.csv'));


    tmax = network.tmax*1e-3;

    % Number of E and I synapses
    N = size(locations,1);
    N_ex_syn = floor(network.parameters.eiFraction*N);
    N_in_syn = N-N_ex_syn;

    % Initialize
    ei = [zeros(1,N_ex_syn),ones(1,N_in_syn)];
    idsE = -ones(ceil(2*N_ex_syn*tmax),1);
    tsE = -ones(ceil(2*N_ex_syn*tmax),1);

    % Network topology
    dist_metric = @(x,y)vecnorm(x-y,2,2);
    elevation_E = locations(1:N_ex_syn,1);
    azimuth_E = locations(1:N_ex_syn,2);
    elevation_I = locations(N_ex_syn+1:N,1);
    azimuth_I = locations(N_ex_syn+1:N,2);

    % Look up table for node indices
    for i = 1:N_ex_syn
        CLUT{i} = find(C(:,1)==i);
    end

    % Simulation parameters
    dt = 4e-3;
    t = 0:dt:tmax;
    tN = length(t);

    % Random external input
    nTrans = poissrnd(network.parameters.eCellParams.firingRate*dt*N_ex_syn);
    idsE(1:nTrans) = randperm(N_ex_syn,nTrans);
    post = idsE(1:nTrans);
    k = nTrans;
    tsE(1:nTrans) = t(1)*ones(nTrans,1);
    count = nTrans;

    exN = poissrnd(network.parameters.eCellParams.firingRate*dt*N_ex_syn*(1-network.branchNo),tN,1);

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
        idsE(count+1:count+exN(i)) = randperm(N_ex_syn,exN(i));
        tsE(count+1:count+exN(i)) = t(i)+dt*rand(exN(i),1);
        count = count+exN(i);
        nTrans = nTrans + exN(i);
    end
    tsE(count:end) = [];
    idsE(count:end) = [];

    % Look up table for excitatory spike times
    speLUT = cell(N_ex_syn,1);
    for j = 1:N_ex_syn
        speLUT{j} = find(idsE==j);
    end

    % Inhibitory neurons follow nearby excitatory neurons
    idsI = -ones(ceil(10*N_in_syn*tmax),1);
    tsI = -ones(ceil(10*N_in_syn*tmax),1);
    count = 0;

    for i = 1:N_in_syn
        idcsJ = C(find(C(:,1)==(N_ex_syn+i)),2)';
        for j = idcsJ
            % Choose excitatory neuron nearby
            idcs = speLUT{j};

            nTrans = length(idcs);
            idsI(count+1:count+nTrans) = i+N_ex_syn;,

            % Make some spikes random depending on branchNo
            idcs2 = rand(nTrans,1)<network.branchNo;
            t0 = tsE(idcs).*idcs2+tmax*rand(nTrans,1).*(1-idcs2);

            % Add spikes to inhibitory neuron
            tsI(count+1:count+nTrans) = t0+dt*rand(nTrans,1);
            count = count+nTrans;
        end
    end
    tsI(count:end) = [];
    idsI(count:end) = [];

    % Combine
    ids = [idsE;idsI];
    ts = [tsE;tsI]*1e3;
end