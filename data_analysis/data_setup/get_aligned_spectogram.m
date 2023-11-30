elName = {'Fz', 'Cz', 'Pz', 'C3', 'C4', 'CP3', 'CP4'};
saveFolder = 'E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\analyzed\rescaled_time';

dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';
load(fullfile(dataFolder,'data_time_information.mat'));
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\raw\time_series_all_channels.mat')

for j = 2
    for i = 1:14
        X = TimeDomainAligned(:,j,i);
        % X = fillgaps(X,1e4);
        idcs = find(and(Time>=-300,Time<60));
        [freq,time,psd(:,:,i)] = eegfft(Time(idcs),X(idcs),2,1.9);

        idcs = find(and(freq>55,freq<65));
        psd(idcs,:,:) = nan;
        idcs = find(and(freq>115,freq<125));
        psd(idcs,:,:) = nan;
        psd(:,:,i) = 10.^fillgaps(log10(squeeze(psd(:,:,i))),5);

        baseline(:,1,i) = nanmedian(psd(:,time<t0(i),i),2);
    end
    idcs = find(freq<150);
    freq = freq(idcs);
    baseline = baseline(idcs,:,:);
    psd = psd(idcs,:,:);

    pre = 10*log10(psd./baseline);
    save(fullfile(saveFolder,sprintf('electrode%d_%s.mat',j,elName{j})),'freq','time','psd','baseine','pre');
end
