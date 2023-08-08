d1 = load('E:\Research_Projects\004_Propofol\data\simulations\raw\test\EI_network\simulation\simulation_data.mat');
eeg1 = network_simulation_beluga.getEEG(detrend(d1.dipoles),sa,1e3);
[freq1,~,psd1] = eegfft(d1.time*1e-3,detrend(eeg1(:,1)),2,1.9);

d2 = load('E:\Research_Projects\004_Propofol\data\simulations\raw\test\EI_network(m=0)\simulation\simulation_data.mat');
eeg2 = network_simulation_beluga.getEEG(detrend(d2.dipoles),sa,1e3);
[freq2,~,psd2] = eegfft(d2.time*1e-3,detrend(eeg2(:,1)),2,1.9);

d3 = load('E:\Research_Projects\004_Propofol\data\simulations\raw\test\EI_network\simulation_tau\simulation_data.mat');
eeg3 = network_simulation_beluga.getEEG(detrend(d3.dipoles),sa,1e3);
[freq3,~,psd3] = eegfft(d3.time*1e-3,eeg3(:,1),2,1.9);

d4 = load('E:\Research_Projects\004_Propofol\data\simulations\raw\test\EI_network\simulation_gamE\simulation_data.mat');
eeg4 = network_simulation_beluga.getEEG(detrend(d4.dipoles),sa,1e3);
[freq4,~,psd4] = eegfft(d4.time*1e-3,eeg4(:,1),2,1.9);

clrs = clrsPT.sequential(10); clrs = clrs(5:end,:);
figureNB(14,8);
subplot(2,3,3+1);
    plot(freq1,mean(psd1,2),'LineWidth',1,'color',clrs(end,:));
    hold on;
    plot(freq2,mean(psd2,2),'LineWidth',1,'color',clrs(1,:));
    xlim([0.5,150])
    ylim([1e-18,1e-14])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xticks([1,10,100]);
    gcaformat
subplot(2,3,3+2);
    plot(freq1,mean(psd1,2),'color',[1,0.6,0.6],'LineWidth',1);
    hold on;
    plot(freq3,mean(psd3,2),'color',[1,0,0],'LineWidth',1);
    xlim([0.5,150])
    ylim([1e-18,1e-14])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xticks([1,10,100]);
    gcaformat
subplot(2,3,3+3);
    plot(freq1,mean(psd1,2),'color',[1,0.6,0.6],'LineWidth',1);
    hold on;
    plot(freq3,mean(psd4,2),'color',[1,0,0],'LineWidth',1);
    xlim([0.5,150])
    ylim([1e-18,1e-14])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xticks([1,10,100]);
    gcaformat



% figureNB(14,4);
subplot(2,2,1);
    raster(ids,ts,gcf);
    set(get(gca,'yaxis'),'visible','on')
    yticks([0,30e3])
    ylim([0,30e3])
    yticklabels({'1','30,000'})
    ylabel('Presynaptic neurons');
    xlabel('Time (ms)');
    gcaformat;
subplot(2,2,2);
    plot(d1.time,eeg0(:,1),'r');
    hold on;
    plot(d1.time,eeg0(:,2),'k');
    xlabel('Time (ms)')
    ylabel(['Single-neuron EEG (' char(956) 'V)'])
    gcaformat


idcs = sa.cortex2K.in_from_cortex75K;
d3 = load('E:\Research_Projects\004_Propofol\data\simulations\raw\test\EI_network\simulation_tau\simulation_data.mat');
eeg3 = network_simulation_beluga.getEEG(d3.dipoles,sa,1e3);
[freq3,~,psd3] = eegfft(d3.time*1e-3,eeg3(:,1),2,1.9);

dp1 = detrend(resample(d1.dipoles(:,:,1),1e3,16e3),'constant');
dp2 = detrend(resample(d2.dipoles(:,:,1),1e3,16e3),'constant');
dp3 = detrend(resample(d3.dipoles(:,:,1),1e3,16e3),'constant');
t = time(1:16:end);
eeg1 = zeros(size(dp1,1),length(idcs));
eeg2 = zeros(size(dp2,1),length(idcs));
eeg3 = zeros(size(dp3,1),length(idcs));
for i = 1:length(idcs)
    eeg1(:,i) = network_simulation_beluga.getEEG(dp1,sa,idcs(i));
    eeg2(:,i) = network_simulation_beluga.getEEG(dp2,sa,idcs(i));
    eeg3(:,i) = network_simulation_beluga.getEEG(dp3,sa,idcs(i));
end

[psd1,f1] = pmtm(eeg1,2,[],1e3);
[psd2,f3] = pmtm(eeg2,2,[],1e3);
[psd3,f3] = pmtm(eeg3,2,[],1e3);

figureNB(4.4,4);
    plot(f1,mean(psd1,2));
    hold on;
    plot(f3,mean(psd3,2));
    xlim([0.5,150])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    gcaformat

figureNB(4.4,4);
    plot(freq2,mean(psd2,2),'LineWidth',1,'color',clrs(1,:));
    hold on;
    plot(freq1,mean(psd0,2),'LineWidth',1,'color',clrs(end-1,:));
    xlim([0.5,150])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    gcaformat