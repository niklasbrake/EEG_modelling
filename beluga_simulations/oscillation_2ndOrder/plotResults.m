load('E:\Research_Projects\004_Propofol\data\simulations\raw\osc_2ndOrder\simulation\simulation_data.mat')
prop = load('E:\Research_Projects\004_Propofol\data\simulations\raw\osc_2ndOrder\simulation_prop\simulation_data.mat');

[freq,~,psd] = eegfft(time(end/2:end)*1e-3,dipoles(end/2:end,3),4,2);
[freq_prop,~,psd_prop] = eegfft(prop.time(end/2:end)*1e-3,prop.dipoles(end/2:end,3),4,2);

[ids,ts,ei] = network.getprenetwork(network.spikingFile);
ids(ts<5e3) = [];
ts(ts<5e3) = [];
k = max(ids)/2;
idcs = find(ids<k);
h = histcounts(ts(idcs),'BinWidth',1);
[freq_h,~,psd_h] = eegfft((1:length(h))'*1e-3,h(:),4,2);

idcs = find(ids>=k);
h = histcounts(ts(idcs),'BinWidth',1);
[freq_h2,~,psd_h2] = eegfft((1:length(h))'*1e-3,h(:),4,2);

figureNB;
subplot(2,1,1);
    R = raster(ids(ids<k),ts(ids<k),gcf);
    R.Color = 'blue';
    hold on;
    R = raster(ids(ids>=k),ts(ids>=k),gcf);
    R.Color = 'green';
    ylim([1,30e3]);
    title('Network activity ~ 2nd order linear system')
subplot(2,2,3);
    plot(freq_h,mean(psd_h,2),'LineWidth',1,'color','b')
    % hold on;
    % plot(freq_h2,mean(psd_h2,2),'LineWidth',1,'color','g')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.5,50])
    xticks([0.5,5,50]);
    xlabel('Frequency (Hz)');
    yticks([]);
    title('Network activity')
    ylabel('log PSD');
subplot(2,2,4);
    plot(freq,mean(psd,2),'color','k','LineWidth',1);
    hold on;
    plot(freq_prop,mean(psd_prop,2),'color','r','LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.5,50])
    xticks([0.5,5,50]);
    xlabel('Frequency (Hz)');
    ylabel('log PSD');
    yticks([]);
    title('EEG')
    L = legend({'\tau = 10 ms','\tau = 30 ms'},'Box','off');
    L.ItemTokenSize = [5,5];
gcaformat(gcf);