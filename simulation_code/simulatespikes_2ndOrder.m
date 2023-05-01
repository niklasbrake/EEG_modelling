function [ids,ts,ei,lam] = simulatespikes_2ndOrder(N,network,oscillation_freq)

dt = 1e-5;
T = 0:dt:network.tmax*1e-3;
wn = 2*pi*oscillation_freq;

a = 0.25;

Y = zeros(2,length(T));
for i = 1:length(T)-1
    Y(1,i+1) = Y(1,i) + dt*Y(2,i);
    Y(2,i+1) = Y(2,i) + dt*(wn^2*randn - 2*a*wn*Y(2,i) - wn^2*Y(1,i));
end
lam = resample(Y(2,:),1e3,round(1/dt));
t = (0:1e-3:10)*1e3;

amplitude = 0.9;
lam = (1-amplitude)+2*amplitude*(lam-min(lam))/(max(lam)-min(lam));


X1 = cumsum(lam)/sum(lam);
    M = floor(N/2);
    mE = network.parameters.eiFraction*M;
    lamE = poissrnd(network.parameters.eCellParams.firingRate*network.tmax*1e-3*mE);
    mI = M-mE;
    lamI = poissrnd(network.parameters.iCellParams.firingRate*network.tmax*1e-3*mI);
    tsE1 = interp1(X1,t,rand(lamE,1),'next','extrap');
    idsE1 = randi(mE,lamE,1);
    tsI1 = interp1(X1,t,rand(lamI,1),'next','extrap');
    idsI1 = randi(mI,lamI,1)+mE;
    ei = [zeros(mE,1);ones(mI,1)];
    elevation1 = rand(M,1)/2+0.5;
    azimuth1 = rand(M,1);

X2 = cumsum(2-lam)/sum(2-lam);
    mE = network.parameters.eiFraction*(N-M);
    lamE = poissrnd(network.parameters.eCellParams.firingRate*network.tmax*1e-3*mE);
    mI = N-M-mE;
    lamI = poissrnd(network.parameters.iCellParams.firingRate*network.tmax*1e-3*mI);
    tsE2 = interp1(X2,t,rand(lamE,1),'next','extrap');
    idsE2 = randi(mE,lamE,1)+M;
    tsI2 = interp1(X2,t,rand(lamI,1),'next','extrap');
    idsI2 = randi(mI,lamI,1)+mE+M;
    ei = [ei;zeros(mE,1);ones(mI,1)];
    elevation2 = rand(N-M,1)/2;
    azimuth2 = rand(N-M,1);


ids = [idsE1;idsI1;idsE2;idsI2];
ts = [tsE1;tsI1;tsE2;tsI2];
elevation = [elevation1;elevation2];
azimuth = [azimuth1;azimuth2];


return;
csvwrite(fullfile(network.preNetwork,'locations.csv'),[azimuth(:),elevation(:)]);
network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,network.spikingFile)
csvwrite(fullfile(network.preNetwork,'multisynapse_IDs.csv'),-ones(N,1));