d1 = load('E:\Research_Projects\004_Propofol\data\simulations\raw\hybrid\simulation\simulation_data.mat');

d2 = load('E:\Research_Projects\004_Propofol\data\simulations\raw\hybrid\simulation_tau\simulation_data.mat');

idcs = sa.cortex2K.in_from_cortex75K;

dp1 = detrend(resample(d1.dipoles(:,:,1),1e3,16e3),'constant');
dp2 = detrend(resample(d2.dipoles(:,:,1),1e3,16e3),'constant');

eeg1 = zeros(size(dp1,1),length(idcs));
eeg2 = zeros(size(dp2,1),length(idcs));
for i = 1:length(idcs)
    eeg1(:,i) = network_simulation_beluga.getEEG(dp1,sa,idcs(i));
    eeg2(:,i) = network_simulation_beluga.getEEG(dp2,sa,idcs(i));
end

[psd1,freq1] = pmtm(eeg1,2,[],1e3);
[psd2,freq2] = pmtm(eeg2,2,[],1e3);


figureNB;
    plot(freq1,mean(psd1,2),'color',[1,0.6,0.6],'LineWidth',1);
    hold on;
    plot(freq2,mean(psd2,2),'color',[1,0,0],'LineWidth',1);
    xlim([0.5,150])
    % ylim([1e-18,1e-14])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xticks([1,10,100]);
    gcaformat


time = d1.time(1:16:end)*1e-3;
waitbar(0);
for i = 1:length(idcs)
    waitbar(i/length(idcs));
    [f1,~,P1_temp] = eegfft(time,eeg1(:,i),2,0);
    P1(:,i) = mean(P1_temp,2);
    [f2,~,P2_temp] = eegfft(time,eeg2(:,i),2,0);
    P2(:,i) = mean(P2_temp,2);
end

figureNB;
    plot(f1,mean(P1,2),'color',[1,0.6,0.6],'LineWidth',1);
    hold on;
    plot(f2,mean(P2,2),'color',[1,0,0],'LineWidth',1);
    xlim([0.5,150])
    % ylim([1e-18,1e-14])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xticks([1,10,100]);
    gcaformat
