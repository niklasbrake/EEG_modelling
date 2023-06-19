passive = load('E:\Research_Projects\005_Aperiodic_EEG\data\simulations\active\LFPy_passive\simulation_data.mat')
active = load('E:\Research_Projects\005_Aperiodic_EEG\data\simulations\active\LFPy_active\simulation_data.mat')

idcs = find(passive.time>100);

passive.time = passive.time(idcs);
passive.V = passive.V(idcs,:);
passive.dipoles = passive.dipoles(idcs,:,:);

active.time = active.time(idcs);
active.V = active.V(idcs,:);
active.dipoles = active.dipoles(idcs,:,:);

[sa,X] = network_simulation_beluga.getHeadModel;
location = randi(size(sa.cortex75K.vc,1)); % Random location

passive.eeg = network_simulation_beluga.getEEG(passive.dipoles,sa,location);
[passive.psd,f] = pmtm(detrend(passive.eeg),2,[],16e3);
passive.psd = passive.psd(2:end,:); f = f(2:end);

active.eeg = network_simulation_beluga.getEEG(active.dipoles,sa,location);
[active.psd,f] = pmtm(detrend(active.eeg),2,[],16e3);
active.psd = active.psd(2:end,:); f = f(2:end);


figureNB;
    plotwitherror(f,passive.psd,'Q');
    hold on;
    plotwitherror(f,active.psd,'Q');
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlabel('Frequency (Hz)');
    xlim([0.5,512]);
    gcaformat;


AP_waveform = [];
for i = 1:10
    [x,y] = findpeaks(active.eeg(:,i),'MinPeakHeight',1e-6);
    for j = 1:length(x)
        waveform_idcs = y(j)-499:y(j)+500;
        waveform_idcs(waveform_idcs<1) = nan;
        waveform_idcs(waveform_idcs>size(active.eeg,1)) = nan;
        subIdcs = find(~isnan(waveform_idcs));
        tempform = nan(1e3,1);
        tempform(subIdcs) = active.eeg(waveform_idcs(subIdcs),i);
        AP_waveform = [AP_waveform,tempform];
    end
end

AP_waveform = nanmean(AP_waveform-nanmedian(AP_waveform),2);
AP_t = [-499:500]/16;
AP_waveform = interp1(AP_t,AP_waveform,AP_t(1):0.1:AP_t(end));
AP_t = AP_t(1):0.1:AP_t(end);


% Initialize network
network = network_simulation_beluga(fullfile('C\Users\brake\Documents\temp','temp_network'));

% Initialize post network
nPostNeurons = 1;
network = network.initialize_postsynaptic_network(nPostNeurons);

% Presyanptic network parameters
nPreNeurons = 30e3;
network.tmax = 2e3; % 2 seconds
network.branchNo = 0.98;
[ids,ts,~,B,t] = network.simulatespikes(0.98);

tvec = 0:0.1:network.tmax;
M = nPreNeurons;
X = zeros(M,length(tvec));
waitbar(0);
for j = 1:M
    waitbar(j/M);
    idcs = interp1(tvec,1:length(tvec),ts(ids==j),'nearest','extrap');
    X(j,idcs) = 1;
    X(j,:) = filter(AP_waveform,1,X(j,:));
end


folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\dyad_network_criticality\m=0.98';
psd = zeros(32769,1);
for j = 1:10
    data1 = csvread(fullfile(folder,['run' int2str(j-1)],'s=0\LFPy\cell00001.csv'),1,0);
    data2 = csvread(fullfile(folder,['run' int2str(j-1)],'s=0\LFPy\cell00002.csv'),1,0);
    eeg1 = network_simulation_beluga.getEEG(data1(2:end-1600,2:4),sa,location);
    eeg2 = network_simulation_beluga.getEEG(data2(2:end-1600,2:4),sa,location);
    psd = psd + pmtm(detrend(eeg1),2,[],16e3) + pmtm(detrend(eeg2),2,[],16e3);
end
psd = psd/20;
[~,f] = pmtm(detrend(eeg2),2,[],16e3);

rho = 0.6;
x0 = samplePSD(f,psd,16e3,1);
x1 = samplePSD(f,psd,16e3,100);
x = rho*x0 + sqrt(1-rho^2)*x1;
