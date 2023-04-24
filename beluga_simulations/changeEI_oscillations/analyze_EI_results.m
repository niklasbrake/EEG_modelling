[sa,X] = network_simulation_beluga.getHeadModel;

folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\EI_oscillations';

F = dir(folder);
goodIdcs = find(and(~strcmp({F(:).name},'..'),~strcmp({F(:).name},'.')));
F = F(goodIdcs);
[~,I] = sort(cellfun(@(x)str2num(x),{F(:).name}));
F = F(I);

locations = randi(size(sa.cortex75K.vc,1),20,1); % Random location

clrs = clrsPT.sequential(length(F)+4);
clrs = clrs(5:end,:);


figureNB;
for j = 1:length(F)
    folder2 = fullfile(folder,F(j).name);
    F2 = dir(folder2);
    F2 = F2(3:end);
    dp = [];
    for i = 1:length(F2)
        try
            data = load(fullfile(folder2,F2(i).name,'LFPy','simulation_data.mat'));
            dipoles = resample(data.dipoles(2:end-1600,:),2e3,16e3);
        catch
            dipoles = nan(20000,3);
        end
        dp(:,:,i) = dipoles;
    end
    % psd{j} = mypmtm(detrend(squeeze(dp(:,2,:)),'constant'),2e3,10);
    psd{j} = zeros(10e3,length(F2));
    for k = 1:length(locations)
        eeg = network_simulation_beluga.getEEG(dp,sa,locations(k));
        psd{j} = psd{j} + mypmtm(detrend(eeg,'constant'),2e3,10)/length(locations);
    end
    plotwitherror(0.1:0.1:1e3,psd{j},'CI','color',clrs(j,:),'LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.1,100]);
    hold on;
    drawnow;
end

f = 0.1:0.1:1e3;

figureNB
    plotwitherror(f,psd{end},'CI','color','k','LineWidth',1);
    plotwitherror(f,psd{end-2},'CI','color','r','LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.1,100]);
    hold on;