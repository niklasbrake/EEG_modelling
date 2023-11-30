dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';

load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\raw\time_series_all_channels.mat')
load(fullfile(dataFolder,'data_time_information.mat'));
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

elName = {'Fz', 'Cz', 'Pz', 'C3', 'C4', 'CP3', 'CP4'};
saveFolder = 'E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\analyzed\rescaled_time';

for j = 1:7
    tCommon = linspace(-1.5,0.5,200);
    for i = 1:14
        X = TimeDomainAligned(:,j,i);
        % X = fillgaps(X,1e4);
        tRescaled = -Time./t0(i);
        idcs = find(and(tRescaled>=-1.6,tRescaled<=0.6));
        [freq,time,psd] = eegfft(Time(idcs),X(idcs),2,1.9);
        tRescaled = -time./t0(i);
        idcs = interp1(tCommon,1:length(tCommon),tRescaled,'nearest');
        psd_new(:,:,i) = splitapply(@(x)nanmean(x,2),psd(1:299,:),idcs);
    end
    psd = psd_new;
    freq = freq(1:299);
    idcs = find(and(freq>55,freq<65));
    psd(idcs,:,:) = nan;
    idcs = find(and(freq>115,freq<125));
    psd(idcs,:,:) = nan;
    for i = 1:200
        psd(:,i,:) = 10.^fillgaps(log10(squeeze(psd(:,i,:))),5);
    end
    tRescaled = tCommon;
    save(fullfile(saveFolder,sprintf('electrode%d_%s.mat',j,elName{j})),'freq','tRescaled','psd');
end