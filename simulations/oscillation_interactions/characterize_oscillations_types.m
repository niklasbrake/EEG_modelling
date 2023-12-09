for i = 1:10
    [eeg0(:,i),eeg1(:,i),eeg2(:,i),eeg3(:,i),eeg4(:,i)] = oscillation_types;
end
[psd0,f] = pmtm(detrend(eeg0,'constant'),2,[],200);
[psd1,f] = pmtm(detrend(eeg1,'constant'),2,[],200);
[psd3,f] = pmtm(detrend(eeg3,'constant'),2,[],200);
[psd4,f] = pmtm(detrend(eeg4,'constant'),2,[],1000);

t  = (0:(length(eeg1)-1))*5;
csvwrite('rhythm1.csv',[t(:),eeg1(:,1)]);
if(any(isnan(eeg2(:,1))))
    warning('ising model output not saved.');
    psd2 = psd0*nan;
else
    csvwrite('rhythm2.csv',[t(:),eeg2(:,1)]);
    [psd2,f] = pmtm(detrend(eeg2,'constant'),2,[],200);
end
csvwrite('rhythm3.csv',[t(:),eeg3(:,1)]);
t4  = (0:(length(eeg4)-1));
csvwrite('rhythm4.csv',[t4(:),eeg4(:,1)]);


psd = {psd0,psd1,psd2,psd3,psd4};
% save('E:\Research_Projects\004_Propofol\data\simulations\raw\peak_trend_sensitivity\rhythm_spectra.mat','f','psd');

clrs = clrsPT.lines(4);

t  = (0:(length(eeg1)-1))*5;
figureNB;
subplot(5,1,1);
    plot(t,eeg0(:,1),'color','k','LineWidth',1);
    % axis off;
    ylim([0.5,1.5]);
    xlim([1e4,2e4]);
subplot(5,1,1+1);
    plot(t,eeg4(:,1),'color',clrs(4,:),'LineWidth',1);
    % axis off;
    ylim([0.5,1.5]);
    xlim([1e4,2e4]);
    % xlim([0,1e4]);
    % xlim([0,1e4]);
subplot(5,1,2+1);
    plot(t,eeg3(:,1),'color',clrs(3,:),'LineWidth',1);
    % axis off;
    ylim([0.5,1.5]);
    xlim([1e4,2e4]);
subplot(5,1,3+1);
    plot(t,eeg2(:,1),'color',clrs(2,:),'LineWidth',1);
    % axis off;
    ylim([0.5,1.5]);
    xlim([1e4,2e4]);
    % xlim([0,1e4]);
subplot(5,1,4+1);
    plot(t,eeg1(:,1),'color',clrs(1,:),'LineWidth',1);
    % axis off;
    ylim([0.5,1.5]);
    xlim([1e4,2e4]);
    % xlim([0,1e4]);

figureNB;
    plot(f,mean(psd0,2),'color','k','LineWidth',1);
    hold on;
    plot(f,mean(psd1,2),'color',clrs(1,:),'LineWidth',1);
    plot(f,mean(psd2,2),'color',clrs(2,:),'LineWidth',1);
    plot(f,mean(psd3,2),'color',clrs(3,:),'LineWidth',1);
    plot(f,mean(psd4,2),'color',clrs(4,:),'LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.1,100])
    xlabel('Frequency (Hz)')
    ylabel('PSD')

function [eeg0,eeg1,eeg2,eeg3,eeg4] = oscillation_types
    t = 1:1e4;

    % Sin waves
    L1 = sin(2*pi*t*10*1e-3);

    % Ising model of oscillations
    folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\test';
    network = network_simulation_beluga(folder);
    fun = fullfile(network.functionFolder,'ising.exe');
    if(~exist(fun))
        warning('Ising model needs to be downloaded from https://doi.org/10.5281/zenodo.7426504, moved to model folder, and compiled to ising.exe');
        eeg2 = nan(size(t));
    else
        beta = 0.98;
        c = 0.005;
        arg = [fun ' ' num2str(beta) ' ' num2str(c) ' > ' fullfile(network.preNetwork,'spikes1.txt')];
        system(arg);
        X1 = dlmread(fullfile(network.preNetwork,'spikes1.txt'));
        t = 1:length(X1);
        L2 = X1(:,2);
        eeg2 = max(1+0.1*L2/std(L2),0);% + randn(size(L2))*1e-2;
    end

    % 2nd order damped oscillation
    dt = 1e-5;
    T = 0:dt:t(end)*1e-3;
    wn = 2*pi*10;
    a = 0.2;
    Y = zeros(2,length(T));
    for i = 1:length(T)-1
        Y(1,i+1) = Y(1,i) + dt*Y(2,i);
        Y(2,i+1) = Y(2,i) + dt*(wn^2*randn - 2*a*wn*Y(2,i) - wn^2*Y(1,i));
    end
    L3 = interp1(T,Y(2,:),t*1e-3);


    [B,t2] = subcritical(0.98);
    L4 = interp1(t2,B,t*1e-3,'linear','extrap')-mean(B);


    eeg0 = max(1+0.1*randn(size(L1)),0);% + randn(size(L1))*1e-2;
    eeg1 = max(1+0.1*L1/std(L1),0);% + randn(size(L1))*1e-2;
    eeg3 = max(1+0.1*L3/std(L3),0);% + randn(size(L3))*1e-2;
    eeg4 = max(1+0.1*L4/std(L4),0);% + randn(size(L4))*1e-2;
end

function [B,t] = subcritical(m)
    tmax = 50;
    M = 1e5;
    dt = 4e-3;
    t = 0:dt:tmax;
    N = length(t);

    k = 4;
    h = 0.75*dt*M*(1-m);

    % Generate mean firing rate using critical branching process
    B = zeros(N,1);
    exN = poissrnd(h,N,1);
    B(1) = exN(1);
    for i = 1:N-1
        if(B(i))
            count = binornd(B(i)*k,m/k);
        else
            count = 0;
        end
        B(i+1) = count+exN(i);
    end
    B(t<2+dt) = [];
    t(t<2+dt) = []; t = t-2;
end