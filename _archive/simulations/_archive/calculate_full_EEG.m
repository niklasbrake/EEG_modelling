
id = 1;
file = ['E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\86955_simulation' sprintf('%02d',id) '.mat'];
while exist(file)
    data = load(file);
    Qbaseline(id).X = data.simulatenous.Qx;
    Qbaseline(id).Y = data.simulatenous.Qy;
    Qbaseline(id).Z = data.simulatenous.Qz;
    id = id+1;
    file = ['E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\86955_simulation' sprintf('%02d',id) '.mat'];
end


id = 1;
file = ['E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\86955_propofol_simulation' sprintf('%02d',id) '.mat'];
while exist(file)
    data2 = load(file);
    Qpropofol(id).X = data.simulatenous.Qx;
    Qpropofol(id).Y = data.simulatenous.Qy;
    Qpropofol(id).Z = data.simulatenous.Qz;
    id = id+1;
    file = ['E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\86955_propofol_simulation' sprintf('%02d',id) '.mat'];
end

neurons = 16e9;


dt = mean(diff(data.t));

M = length(Qbaseline);
N = 1000;
P = [];
for j = 1:30
    W = zeros(length(Qbaseline(1).X),N);
    for i = 1:N
        id = randi(M);
        W(:,i) = getEEG(Qbaseline(id),sa);
    end

    eeg = mean(W,2)*neurons;
    [P(:,j),F] = pmtm(eeg,2,[],1e3/dt);
    j
end

M = length(Qpropofol);
N = 1000;
P2 = [];
for j = 1:30
    W = zeros(length(Qpropofol(1).X),N);
    for i = 1:N
        id = randi(M);
        W(:,i) = getEEG(Qpropofol(id),sa);
    end

    eeg = mean(W,2)*neurons;
    [P2(:,j),F] = pmtm(eeg,2,[],1e3/dt);
    j
end

figureNB;
plot(F,mean(P,2),'LineWidth',1,'color','k'); hold on;
plot(F,mean(P2,2),'LineWidth',1,'color','r');
set(gca,'xscale','log')
set(gca,'yscale','log')
ylim([1e-4,1e4])
xlim([0.5,300])