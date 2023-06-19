load('E:\Research_Projects\005_Aperiodic_EEG\data\simulations\bAP_unitary_response\unitary_waveform')
AP_t = [-499:500]/16;
AP_waveform = interp1(AP_t,mean(unitary_waveform,2),AP_t(1):0.1:AP_t(end));
AP_t = AP_t(1):0.1:AP_t(end);
%%%%%%%%%%%%%
% Initialize network
network = network_simulation_beluga(fullfile('C\Users\brake\Documents\temp','temp_network'));

% Initialize post network
nPostNeurons = 1;
network = network.initialize_postsynaptic_network(nPostNeurons);

% Presyanptic network parameters
nPreNeurons = 30e3;
network.tmax = 10e3; % 2 seconds


[ids,ts,~,B,t] = network.simulatespikes(m);

tvec = 0:0.1:network.tmax;
X1 = zeros(M,length(tvec));
X2 = zeros(M,length(tvec));
waitbar(0);
for j = 1:M
    waitbar(j/M);
    t0 = ts(ids==j);
    idcs = interp1(tvec,1:length(tvec),t0,'nearest','extrap');
    X1(j,idcs) = 1;
    X1(j,:) = filter(AP_waveform,1,X1(j,:));

    t1 = rand(size(t0))*network.tmax;
    idcs = interp1(tvec,1:length(tvec),t1,'nearest','extrap');
    X2(j,idcs) = 1;
    X2(j,:) = filter(AP_waveform,1,X2(j,:));
end
figureNB(12,7.7);
subplot(2,1,1);
    plot(sum(X1),'k')
    xlim([1,10e4])
    ylim([-1,11]*1e-6);
subplot(2,1,2);
    raster(ids,ts,gcf);


figureNB;

mValues = linspace(0.86,0.99,5);
M = 500;
N = (1:M)';
B2 = [N,ones(M,1)];
B3 = [N.*(N-1)];
V1 = zeros(M,length(mValues));
V2 = zeros(M,length(mValues));
for im = 1:length(mValues)
    m = mValues(im);
    [ids,ts,~,B,t] = network.simulatespikes(m);

    tvec = 0:0.1:network.tmax;
    X1 = zeros(M,length(tvec));
    X2 = zeros(M,length(tvec));
    waitbar(0);
    for j = 1:M
        waitbar(j/M);
        t0 = ts(ids==j);
        idcs = interp1(tvec,1:length(tvec),t0,'nearest','extrap');
        X1(j,idcs) = 1;
        X1(j,:) = filter(AP_waveform,1,X1(j,:));

        t1 = rand(size(t0))*network.tmax;
        idcs = interp1(tvec,1:length(tvec),t1,'nearest','extrap');
        X2(j,idcs) = 1;
        X2(j,:) = filter(AP_waveform,1,X2(j,:));
    end
    waitbar(0);
    for i = 1:M
        waitbar(i/M);
        idcs = randsample(M,i);
        mX1 = sum(X1(idcs,:));
        mX2 = sum(X2(idcs,:));
        V1(i,im) = var(mX1);
        V2(i,im) = var(mX2);
        drawnow;
    end
    
    S0 = N\V2(:,im)
    nullV = S0*N;
    K = B3\(V1(:,im)-nullV);
    rho(im) = K/S0
    subplot(1,2,1);
        plot(N,V2(:,im),'color','k');
        hold on;
        plot(N,nullV,'--b');
        xlim([0,500]);
        xlabel('N');
        ylabel(['\sigma^2 (' char(956) 'V^2)'])
    subplot(1,2,2);
        plot(N,V1(:,im),'color','k');
        hold on;
        plot(N,(N+N.*(N-1)*rho(im))*S0,'--b');
        xlim([0,500]);
        xlabel('N');
        ylabel(['\sigma^2 (' char(956) 'V^2)'])
    drawnow;
end