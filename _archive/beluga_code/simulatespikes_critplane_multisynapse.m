function [ids,ts,ei,C,elevation,azimuth,parents] = simulatespikes_critplane(N,branchNo,tmax)

% Number of E and I synapses
N_ex_syn = floor(0.8*N);
N_in_syn = N-N_ex_syn;

% Multisynapse count (predicted by Markram et al., Cell 2015)
eMulti = 4.5;
iMulti = 13.9;
% iMulti = 4.5;
N_ex_pre = ceil(N_ex_syn/eMulti);
N_in_pre = ceil(N_in_syn/iMulti);

% Initialize
ei = [zeros(1,N_ex_syn),ones(1,N_in_syn)];
idsE = -ones(ceil(2*N_ex_pre*tmax),1);
tsE = -ones(ceil(2*N_ex_pre*tmax),1);
lamE = 0.5;
lamI = lamE*5;

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
nTrans = poissrnd(lamE*dt*N_ex_pre);
idsE(1:nTrans) = randperm(N_ex_pre,nTrans);
post = idsE(1:nTrans);
k = nTrans;
tsE(1:nTrans) = t(1)*ones(nTrans,1);
count = nTrans;

exN = poissrnd(lamE*dt*N_ex_pre*(1-branchNo),tN,1);

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

exN = poissrnd(lamI/lamE*dt*N_in_pre*(1-branchNo),tN,1);
for i = 1:N_in_pre
    D0 = dist_metric([elevation_I(i),azimuth_I(i)],[elevation_E,azimuth_E]);
    [~,I] = sort(D0,'ascend');
    k = poissrnd(lamI/lamE);
    for j = 2:k+1
        idcs = speLUT{I(j)};
        nTrans = length(idcs);
        idsI(count+1:count+nTrans) = i+N_ex_syn;,

        % Make some spikes random depending on branchNo
        idcs2 = rand(nTrans,1)<branchNo;
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
    e = truncate(makedist('Normal','mu',elevation_E(M(j)),'sigma',0.05),0,1);
    a = truncate(makedist('Normal','mu',azimuth_E(M(j)),'sigma',0.05),0,1);

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
