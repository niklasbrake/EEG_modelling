function [ids,ts,ei,x,X1] = simulatespikes_avaETosc(network,beta,c,N)

functionFolder = fileparts(mfilename('fullpath'));
fun = fullfile(functionFolder,'ising.exe');

% Run ising model twice
arg = [fun ' ' num2str(beta) ' ' num2str(c) ' > ' fullfile(network.preNetwork,'spikes1.txt')];
system(arg);
% system([fun ' ' num2str(beta) ' > ' fullfile(network.preNetwork,'spikes2.txt')]);

ei = [];
% Apical network
% Read mean firing rate from adaptive Ising model
    X1 = dlmread(fullfile(network.preNetwork,'spikes1.txt'));
    X2 = -X1;
    t = 0:network.tmax;
    X1 = X1(:,1); X1 = X1(1:network.tmax+1);
    X1 = X1./abs(min(X1))+1+1e-6;
% Randomly draw spike times from this firing rate
    M = floor(N/2);
    mE = network.eiFraction*M;
    lamE = poissrnd(network.eFiringRate*network.tmax*1e-3*mE);
    mI = M-mE;
    lamI = poissrnd(network.iFiringRate*network.tmax*1e-3*mI);
    tsE1 = interp1(cumsum(X1)/sum(X1),t,rand(lamE,1),'next','extrap');
    idsE1 = randi(mE,lamE,1);
    ei = [ei;zeros(mE,1)];
    tsI1 = interp1(cumsum(X1)/sum(X1),t,rand(lamI,1),'next','extrap');
    idsI1 = randi(mI,lamI,1)+mE;
    ei = [ei;ones(mI,1)];
% Random distribution elevation in top half of sphere
    elevation1 = rand(M,1)/2+0.5;
    azimuth1 = rand(M,1);

% Basal network
% Read mean firing rate from adaptive Ising model
    % X2 = dlmread(fullfile(network.preNetwork,'spikes2.txt'));
    X2 = X2(:,1); X2 = X2(1:network.tmax+1);
    X2 = X2./abs(min(X2))+1+1e-6;
% Randomly draw spike times from this firing rate
    mE = network.eiFraction*(N-M);
    lamE = poissrnd(network.eFiringRate*network.tmax*1e-3*mE);
    mI = N-M-mE;
    lamI = poissrnd(network.iFiringRate*network.tmax*1e-3*mI);
    tsE2 = interp1(cumsum(X2)/sum(X2),t,rand(lamE,1),'next','extrap');
    idsE2 = randi(mE,lamE,1)+M;
    ei = [ei;zeros(mE,1)];
    tsI2 = interp1(cumsum(X2)/sum(X2),t,rand(lamI,1),'next','extrap');
    idsI2 = randi(mI,lamI,1)+mE+M;
    ei = [ei;ones(mI,1)];
% Random distribution elevation in bottom half of sphere
    elevation2 = rand(N-M,1)/2;
    azimuth2 = rand(N-M,1);

ids = [idsE1;idsI1;idsE2;idsI2];
ts = [tsE1;tsI1;tsE2;tsI2];
elevation = [elevation1;elevation2];
azimuth = [azimuth1;azimuth2];
x = [azimuth(:),elevation(:)];