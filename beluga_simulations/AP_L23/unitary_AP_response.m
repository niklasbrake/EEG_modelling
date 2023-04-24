[sa,X] = network_simulation_beluga.getHeadModel;
locations = sa.cortex2K.in_from_cortex75K;
locations = randsample(sa.cortex2K.in_from_cortex75K,20);

passive = load('E:\Research_Projects\005_Aperiodic_EEG\data\simulations\active\LFPy_passive\simulation_data.mat');
active = load('E:\Research_Projects\005_Aperiodic_EEG\data\simulations\active\LFPy_active\simulation_data.mat');

idcs = find(passive.time>100);
K = size(passive.V,2);

passive.time = passive.time(idcs);
passive.V = passive.V(idcs,:);
passive.dipoles = passive.dipoles(idcs,:,:);
passive.eeg = network_simulation_beluga.getEEG(passive.dipoles,sa,48108);

active.time = active.time(idcs);
active.V = active.V(idcs,:);
active.dipoles = active.dipoles(idcs,:,:);
active.eeg = network_simulation_beluga.getEEG(active.dipoles,sa,48108);

for i = 1:K
    [x,y] = findpeaks(active.eeg(:,i)-passive.eeg(:,i),'MinPeakHeight',1e-6);
    t{i} = y;
end

X_AP = zeros(1e3,length(locations));
P_passive = zeros(16385,length(locations));
P_active = zeros(16385,length(locations));
for k = 1:length(locations)
    waitbar(k/length(locations))
    passive.eeg = network_simulation_beluga.getEEG(passive.dipoles,sa,locations(k));
    active.eeg = network_simulation_beluga.getEEG(active.dipoles,sa,locations(k));

    AP_waveform = zeros(1000,1);
    for i = 1:K
        y = t{i};
        for j = 1:length(t{i})
            waveform_idcs = y(j)-499:y(j)+500;
            waveform_idcs(waveform_idcs<1) = nan;
            waveform_idcs(waveform_idcs>size(active.eeg,1)) = nan;
            subIdcs = find(~isnan(waveform_idcs));
            tempform = zeros(1e3,1);
            tempform(subIdcs) = active.eeg(waveform_idcs(subIdcs),i)-passive.eeg(waveform_idcs(subIdcs),i);
            AP_waveform = AP_waveform + tempform;
        end
    end
    X_AP(:,k) = AP_waveform;

    P_passive(:,k) = mean(pmtm(detrend(passive.eeg),2,[],16e3),2);
    P_active(:,k) = mean(pmtm(detrend(active.eeg),2,[],16e3),2);
end
[~,f] = pmtm(detrend(passive.eeg),2,[],16e3);
f = f(2:end);

s = floor(size(passive.eeg,1)/2);
temp = zeros(size(passive.eeg,1),k);
temp(s-499:s+500,:) = (X_AP-median(X_AP))/sum(cellfun(@(x)length(x),t));
[psd0,f0] = pmtm(temp,2,[],16e3);
psd0 = psd0(2:end,:); f0 = f0(2:end);

figureNB(6.5,8);
axes('Position',[0.185,0.71,0.775,0.21]);
    plot(mean(X_AP,2),'color','k','LineWidth',1);
    ylabel(['EEG (' char(956) 'V)'])
    ylim([-15,15]*1e-7)
    set(get(gca,'xaxis'),'visible','off');
    line([800,960],[-1,-1]*1e-6,'color','k','LineWidth',2)
    text(880,-1.2e-6,'10 ms','FontSize',6,'HorizontalAlignment','center','VerticalAlignment','top');
    gcaformat;
axes('Position',[0.185,0.11,0.775,0.515]);
    plotwitherror(f,P_passive(2:end,:),'CI');
    hold on;
    plotwitherror(f,P_active(2:end,:),'CI');
    plotwitherror(f0,psd0,'CI','color','k');
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlabel('Frequency (Hz)');
    xlim([0.5,8e3])
    ylim([1e-25,1e-15])
    gcaformat;



% Compute expected EEG variance of 16 billion passive neurons
dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';
load(fullfile(dataFolder,'anatomy_cortical_pairwise_distance_distribution.mat'));
signed_area = A;
total_area = B;
N = 16e9;
dMids = 0.5*(rValues(2:end)+rValues(1:end-1));
nrnCount = mean(diff(signed_area),2)*200000;
nrnCount(end) = N-sum(nrnCount(1:end-1));
corr_kernel = @(d) exp(-d.^2/4);
rho_bar = sum(corr_kernel(dMids).*nrnCount)/sum(nrnCount');
SIG_N = @(rho) N+N*(N-1)*rho;


figureNB;
subplot(2,1,1);
    plot(f,mean(P_passive(2:end,:),2))
    hold on;
    plot(f0,mean(psd0,2));
    plot(f0,mean(psd0,2)+mean(P_passive(2:end,:),2),'color','k');
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlabel('Frequency (Hz)');
    xlim([0.5,8e3])
    ylim([1e-25,1e-15])
    gcaformat;
subplot(2,1,2);
    syn_scaled = SIG_N(0.14*rho_bar)*mean(P_passive(2:end,:),2);
    ap_scaled = SIG_N(0.008*rho_bar)*mean(psd0,2);
    plot(f,syn_scaled)
    hold on;
    plot(f0,ap_scaled);
    plot(f0,syn_scaled+ap_scaled,'color','k');
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlabel('Frequency (Hz)');
    xlim([0.5,8e3])
    gcaformat;


syn_scaled = SIG_N(0.14*rho_bar)*mean(P_passive(2:end,:),2);
T = linspace(0,0.14,1e3);
for i = 1:length(T)
    ap_scaled = SIG_N(T(i)*rho_bar)*mean(psd0,2);
    v(i) = sum((syn_scaled+ap_scaled)./syn_scaled)*mean(diff(f));
end