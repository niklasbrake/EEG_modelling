% for i = 1:10
    [eeg1,eeg2,eeg3] = oscillation_types;
%     [psd1(:,i),f] = pmtm(detrend(eeg1,'constant'),2,[],1e3);
%     [psd2(:,i),f] = pmtm(detrend(eeg2,'constant'),2,[],1e3);
%     [psd3(:,i),f] = pmtm(detrend(eeg3,'constant'),2,[],1e3);
% end
return;
clrs = clrsPT.lines(3);

figureNB;
subplot(3,1,1);
    plot(t,eeg1,'color',clrs(1,:),'LineWidth',1);
    xlim([5e3,6e3]); axis off;
    ylim([0,3]);
subplot(3,1,2);
    plot(t,eeg2,'color',clrs(2,:),'LineWidth',1);
    xlim([5e3,6e3]); axis off;
    ylim([0,3]);
subplot(3,1,3);
    plot(t,eeg3,'color',clrs(3,:),'LineWidth',1);
    xlim([5e3,6e3]); axis off;
    ylim([0,3]);

figureNB;
    plot(f,mean(psd1,2),'color',clrs(1,:),'LineWidth',1);
    hold on;
    plot(f,mean(psd2,2),'color',clrs(2,:),'LineWidth',1);
    plot(f,mean(psd3,2),'color',clrs(3,:),'LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.1,100])
    xlabel('Frequency (Hz)')
    ylabel('PSD')

function [eeg1,eeg2,eeg3] = oscillation_types
    % Ising model of oscillations
    folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\test';
    network = network_simulation_beluga(folder);
    fun = fullfile(network.functionFolder,'ising.exe');
    beta = 0.98;
    c = 0.005;
    arg = [fun ' ' num2str(beta) ' ' num2str(c) ' > ' fullfile(network.preNetwork,'spikes1.txt')];
    system(arg);
    X1 = dlmread(fullfile(network.preNetwork,'spikes1.txt'));
    t = 1:length(X1);
    L2 = X1(:,2);

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

    % Sin waves
    L1 = sin(2*pi*t*10*1e-3);



    eeg1 = max(1+0.4*L1/std(L1),0) + randn(size(L1))*1e-2;
    eeg2 = max(1+0.4*L2/std(L2),0) + randn(size(L2))*1e-2;
    eeg3 = max(1+0.4*L3/std(L3),0) + randn(size(L3))*1e-2;
end
