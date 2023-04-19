function [ids,ts,ei] = simulatespikes_osc(N,tmax,network,oscillation_freq)
    if(nargin<4)
        oscillation_freq = 2; % Hz
    end
    M1 = floor(N/2);
    M2 = N-M1;
    [ids1,ts1,ei1,elevation1,azimuth1,parents1] = simulatespikes(M1,tmax,0,network.branchNo,oscillation_freq);
    [ids2,ts2,ei2,elevation2,azimuth2,parents2] = simulatespikes(M2,tmax,1,network.branchNo,oscillation_freq);
    ids = [ids1;ids2+M1];
    ts = [ts1;ts2];
    ei = [ei1(:);ei2(:)];
    elevation = [elevation1;elevation2];
    azimuth = [azimuth1;azimuth2];
    parents2(parents2~=-1) = parents2(parents2~=-1)+M1;
    parents = [parents1;parents2];
    csvwrite(fullfile(network.preNetwork,'locations.csv'),[azimuth(:),elevation(:)]);
    csvwrite(fullfile(network.preNetwork,'multisynapse_IDs.csv'),parents(:));
    network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,network.spikingFile)
end
function [ids,ts,ei,elevation,azimuth,parents] = simulatespikes(N,tmax,offset,branchNo,f0)

    tmax = tmax*1e-3;
    eFiringRate = 0.5;
    iFiringRate = 2.5;

    % Number of E and I synapses
    N_ex_syn = floor(0.8*N);
    N_in_syn = N-N_ex_syn;

    % Multisynapse count (predicted by Markram et al., Cell 2015)
    % N_ex_pre = ceil(N_ex_syn/3.6);
    % N_in_pre = ceil(N_in_syn/13.9);
    N_ex_pre = ceil(N_ex_syn/1);
    N_in_pre = ceil(N_in_syn/1);

    % Initialize
    ei = [zeros(1,N_ex_syn),ones(1,N_in_syn)];
    idsE = -ones(ceil(2*N_ex_pre*tmax),1);
    tsE = -ones(ceil(2*N_ex_pre*tmax),1);

    % Network topology
    elevation_E = rand(N_ex_pre,1)/2+0.5*offset;
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
        C(nNeigh*(i-1)+1:nNeigh*i,1) = i*ones(nNeigh,1);
        C(nNeigh*(i-1)+1:nNeigh*i,2) = I;
    end
    C(:,3) = branchNo/nNeigh;

    % Look up table for node indices
    for i = 1:N_ex_pre
        CLUT{i} = find(C(:,1)==i);
    end

    % Simulation parameters
    dt = 4e-3;
    t = 0:dt:tmax;
    tN = length(t);

    % Random external input
    nTrans = poissrnd(eFiringRate*dt*N_ex_pre);
    idsE(1:nTrans) = randperm(N_ex_pre,nTrans);
    post = idsE(1:nTrans);
    k = nTrans;
    tsE(1:nTrans) = t(1)*ones(nTrans,1);
    count = nTrans;

    % osc = 1+0.7*sin(2*pi*t*f0+pi*offset);
    osc = 1+sin(2*pi*t*f0+pi*offset);

    for i = 2:tN-1
        % Get spiking cells at previous time point
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
        lamE = eFiringRate*osc(i);
        exN(i) = poissrnd(lamE*dt*N_ex_pre*(1-branchNo));
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
    elevation_I = rand(N_in_pre,1)/2+0.5*offset;
    azimuth_I = rand(N_in_pre,1);
    idsI = -ones(ceil(10*N_in_pre*tmax),1);
    tsI = -ones(ceil(10*N_in_pre*tmax),1);
    count = 0;

    oscDist = cumsum(osc)/sum(osc);

    for i = 1:N_in_pre
        D0 = dist_metric([elevation_I(i),azimuth_I(i)],[elevation_E,azimuth_E]);
        [~,I] = sort(D0,'ascend');
        k = poissrnd(iFiringRate/eFiringRate);
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
            idcs2 = rand(nTrans,1)<branchNo;
            t0 = tsE(idcs).*idcs2+interp1(oscDist,t,rand(nTrans,1),'next','extrap').*(1-idcs2);

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

    a0 = 0.5*offset;
    a1 = 0.5+0.5*offset;
    % Final inhibitory synapses
    M = randi(N_in_pre,N_in_syn-N_in_pre,1);
    for j = 1:length(M)
        % Get random presyanptic neuron
        idcs = spiLUT{M(j)};

        % Truncated normal distribution centered on "parent" synapse
        e = truncate(makedist('Normal','mu',elevation_I(M(j)),'sigma',0.02),a0,a1);
        a = truncate(makedist('Normal','mu',azimuth_I(M(j)),'sigma',0.02),0,1);

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
end