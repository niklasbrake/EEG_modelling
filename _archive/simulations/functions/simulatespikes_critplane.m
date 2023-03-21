function [ids,ts,ei,C,elevation,azimuth] = simulatespikes_critplane(N,branchNo,tmax)

% Number of E and I synapses
N_ex_syn = floor(0.8*N);
N_in_syn = N-N_ex_syn;

% Initialize
ei = [zeros(1,N_ex_syn),ones(1,N_in_syn)];
idsE = -ones(ceil(2*N_ex_syn*tmax),1);
tsE = -ones(ceil(2*N_ex_syn*tmax),1);
lamE = 1;
lamI = 5;

% Network topology
elevation_E = rand(N_ex_syn,1);
azimuth_E = rand(N_ex_syn,1);
dist_metric = @(x,y)vecnorm(x-y,2,2);
nNeigh = 4;
C = zeros(N_ex_syn*nNeigh,3);
for i = 1:N_ex_syn
    D0 = dist_metric([elevation_E(i),azimuth_E(i)],[elevation_E,azimuth_E]);
    [~,I] = sort(D0,'ascend'); I = I(2:nNeigh+1);
    C(nNeigh*(i-1)+1:nNeigh*i,1) = i*ones(nNeigh,1);
    C(nNeigh*(i-1)+1:nNeigh*i,2) = I;
end
C(:,3) = branchNo/nNeigh;

% Look up table for node indices
for i = 1:N_ex_syn
    CLUT{i} = find(C(:,1)==i);
end

% Simulation parameters
dt = 4e-3;
t = 0:dt:tmax;
tN = length(t);

% Random external input
nTrans = poissrnd(lamE*dt*N_ex_syn);
idsE(1:nTrans) = randperm(N_ex_syn,nTrans);
post = idsE(1:nTrans);
k = nTrans;
tsE(1:nTrans) = t(1)*ones(nTrans,1);
count = nTrans;

exN = poissrnd(lamE*dt*N_ex_syn*(1-branchNo),tN,1);

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

% Look up table for spike times
spLUT = cell(N_ex_syn,1);
for j = 1:N_ex_syn
    spLUT{j} = find(idsE==j);
end

% Inhibitory neurons follow nearby excitatory neurons
elevation_I = rand(N_in_syn,1);
azimuth_I = rand(N_in_syn,1);
idsI = -ones(ceil(10*N_in_syn*tmax),1);
tsI = -ones(ceil(10*N_in_syn*tmax),1);
count = 0;
for i = 1:N_in_syn
    D0 = dist_metric([elevation_I(i),azimuth_I(i)],[elevation_E,azimuth_E]);
    [~,I] = sort(D0,'ascend');
    for j = 2:6
        idcs = spLUT{I(j)};
        nTrans = length(idcs);
        idsI(count+1:count+nTrans) = i+N_ex_syn;
        tsI(count+1:count+nTrans) = tsE(idcs)+dt*rand(nTrans,1);
        count = count+nTrans;
    end
end
tsI(count:end) = [];
idsI(count:end) = [];


% Combine
elevation = [elevation_E;elevation_I];
azimuth = [azimuth_E;azimuth_I];
ids = [idsE;idsI];
ts = [tsE;tsI]*1e3;
