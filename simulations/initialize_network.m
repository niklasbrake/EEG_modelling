function [ids,ts,m0] = initialize_network(masterPath,branchNo)
% masterPath = '/lustre04/scratch/nbrake/data/simulations/raw/correlation_test';
addpath('/home/nbrake/aperiodic_EEG_modelling/simulations/functions');

tmax0 = 5e3;
tmax = 2e3;
nrnCount = 2;

network = network_simulation_beluga(masterPath);
network = network.initialize_postsynaptic_network(nrnCount,ones(nrnCount,1));
network.tmax = tmax;
% network = network.setsynapsecount(1e3);
N = network.getsynapsecount;

% Plane
elevation = rand(N,1);
azimuth = rand(N,1);
x = elevation;
y = azimuth;
dist_metric = @(x,y)vecnorm(x-y,2,2);

dt = 4e-3;
nNeigh = 4;
C = zeros(N*nNeigh,3);
for i = 1:N
    D0 = dist_metric([elevation(i),azimuth(i)],[elevation,azimuth]);
    D = exp(-D0.^2/1);
    D(i) = min(D);
    idcs = find(D>1e-9);
    I = [];
    if(length(idcs)<nNeigh)
        [~,I] = sort(D0,'ascend'); I = I(2:nNeigh+1);
    else
        for j = 1:nNeigh
            j0 = interp1(cumsum(D(idcs))/sum(D(idcs)),idcs,rand,'next','extrap');
            I(j) = j0;
            idcs = setdiff(idcs,j0);
        end
    end
    [~,I] = sort(D0,'ascend'); I = I(2:nNeigh+1);
    % I = randperm(N,nNeigh);
    C(nNeigh*(i-1)+1:nNeigh*i,1) = i*ones(nNeigh,1);
    C(nNeigh*(i-1)+1:nNeigh*i,2) = I;
    % C(nNeigh*(i-1)+1:nNeigh*i,4) = dt*2*rand(nNeigh,1); % Random transmission speed
end
C(:,3) = rand(size(C,1),1)*2*branchNo/nNeigh;
C(:,3) = branchNo/nNeigh;

[ids,ts,ei,~,m0] = simulatespikes_det(N,branchNo,tmax0*1e-3,C,nNeigh);

file = fullfile(network.preNetwork,'spikeTimesLong.csv');
network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,file)

correlationFile = fullfile(network.preNetwork,'correlations.csv');
network = network.setCorrelationFile(correlationFile);

csvwrite(fullfile(network.preNetwork,'connections.csv'),C);
csvwrite(fullfile(network.preNetwork,'locations.csv'),[elevation(:),azimuth(:)]);

network.save();