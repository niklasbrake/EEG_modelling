d1 = load('E:\Research_Projects\004_Propofol\data\simulations\raw\hybrid\simulation\simulation_data.mat');

d2 = load('E:\Research_Projects\004_Propofol\data\simulations\raw\hybrid\simulation_tau\simulation_data.mat');

d3 = load('E:\Research_Projects\004_Propofol\data\simulations\raw\hybrid(m=0)\simulation\simulation_data.mat');

d4 = load('E:\Research_Projects\004_Propofol\data\simulations\raw\hybrid(m=0)\simulation_tau\simulation_data.mat');

idcs = sa.cortex2K.in_from_cortex75K;

dp1 = detrend(resample(d1.dipoles(:,:,1),1e3,16e3),'constant');
dp2 = detrend(resample(d2.dipoles(:,:,1),1e3,16e3),'constant');
dp3 = detrend(resample(d3.dipoles(:,:,1),1e3,16e3),'constant');
dp4 = detrend(resample(d4.dipoles(:,:,1),1e3,16e3),'constant');


eeg1 = zeros(size(dp1,1),length(idcs));
eeg2 = zeros(size(dp2,1),length(idcs));
eeg3 = zeros(size(dp3,1),length(idcs));
eeg4 = zeros(size(dp4,1),length(idcs));
for i = 1:length(idcs)
    eeg1(:,i) = network_simulation_beluga.getEEG(dp1,sa,idcs(i));
    eeg2(:,i) = network_simulation_beluga.getEEG(dp2,sa,idcs(i));
    eeg3(:,i) = network_simulation_beluga.getEEG(dp3,sa,idcs(i));
    eeg4(:,i) = network_simulation_beluga.getEEG(dp4,sa,idcs(i));
end

time = d1.time(1:16:end)*1e-3;
waitbar(0);
for i = 1:length(idcs)
    waitbar(i/length(idcs));
    [f1,~,P1_temp] = eegfft(time,eeg1(:,i),2,0);
    P1(:,i) = mean(P1_temp,2);
    [f2,~,P2_temp] = eegfft(time,eeg2(:,i),2,0);
    P2(:,i) = mean(P2_temp,2);
    [f3,~,P3_temp] = eegfft(time,eeg3(:,i),2,0);
    P3(:,i) = mean(P3_temp,2);
    [f4,~,P4_temp] = eegfft(time,eeg4(:,i),2,0);
    P4(:,i) = mean(P4_temp,2);
end



load('E:\Research_Projects\004_Propofol\data\simulations\raw\hybrid\model.mat')

clrs = clrsPT.sequential(10); clrs = clrs(5:end,:);
figureNB(14,8);;
subplot(2,2,1);
    [ids,ts,ei] = network.getprenetwork(fullfile(network.preNetwork,'spikeTimes.csv'));
    raster(ids,ts,gcf);
    set(get(gca,'yaxis'),'visible','on')
    yticks([0,14156,30e3])
    ylim([0,30e3])
    yticklabels({'1','14,156','30,000'})
    ylabel('Neurons');
    title('Full network');
    gcaformat;
subplot(2,2,2);
    [ids,ts,ei] = network.getprenetwork(network.spikingFile);
    raster(ids,ts,gcf)
    set(get(gca,'yaxis'),'visible','on')
    yticks([0,5e3])
    ylim([0,length(ei)])
    yticklabels({'1','5,000'})
    title('Neuron #14156');
    ylabel('Presyanptic neurons');
    gcaformat;
subplot(2,3,4);
    plot(f2,mean(P3,2),'color',clrs(1,:),'LineWidth',1);
    hold on;
    plot(f1,mean(P1,2),'color',clrs(end,:),'LineWidth',1);
    xlim([0.5,150])
    % ylim([1e-18,1e-14])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xticks([1,10,100]);
    ylim([1e-])
    gcaformat
subplot(2,3,5);
    plot(f1,mean(P3,2),'color',[1,0.6,0.6],'LineWidth',1);
    hold on;
    plot(f2,mean(P4,2),'color',[1,0,0],'LineWidth',1);
    xlim([0.5,150])
    % ylim([1e-18,1e-14])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xticks([1,10,100]);
    ylim([1e-])
    gcaformat
subplot(2,3,6);
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
    ylim([1e-])
    gcaformat
