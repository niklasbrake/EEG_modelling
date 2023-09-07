dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';
load(fullfile(dataFolder,'data_time_information.mat'));
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\raw\time_series_all_channels.mat')
for i = 1:14
    X = TimeDomainAligned(:,2,i);
    X = fillgaps(X,1e4);
    idcs = find(and(Time>=-300,Time<60));
    [aligned.freq,aligned.time,aligned.psd(:,:,i)] = eegfft(Time(idcs),X(idcs),2,1.9);
    aligned.baseline(:,1,i) = nanmedian(aligned.psd(:,aligned.time<t0(i),i),2);
end
idcs = find(aligned.freq<150);
aligned.freq = aligned.freq(idcs);
aligned.baseline = aligned.baseline(idcs,:,:);
aligned.psd = aligned.psd(idcs,:,:);
aligned.pre = 10*log10(aligned.psd./aligned.baseline);