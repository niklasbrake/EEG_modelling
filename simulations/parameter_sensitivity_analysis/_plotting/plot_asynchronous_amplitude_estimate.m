load(fullfile(dataFolder,'simulations','simulated_spectra_parameter_sensitivity.mat'));

% Plot histogram of total power
figureNB(3,2.8);
    histogram(log10(sum(psd.*mean(diff(f)))),'EdgeColor','none', ...
        'FaceColor',[1,1,1]*0.4,'FaceAlpha',1)
    str = sprintf(['Single-neuron EEG\naverage power (' char(956) 'V)']);
    xlabel(str);
    yticks([]);
    xlim([-16.5,-11.5]);
    xticks([-16,-14,-12]);
    xticklabels({'10^{-16}','10^{-14}','10^{-12}'});
    gcaformat

data = load(fullfile(dataFolder,'EEG_data','electrode2_Cz.mat'));
data.baseline = squeeze(nanmedian(data.psd(:,data.time<-1,:),2));

% Scale untiary spectra by 16 billion and compare to data
figureNB(3,2.8);
    plotwitherror(data.freq,data.baseline,'M','LineWidth',1,'color','k');
    plotwitherror(f,psd*16e9,'Q','color',[1,1,1]*0.4,'LineWidth',1);
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    set(gca,'xscale','log');
    set(gca,'yscale','log')
    xlim([1,100])
    xticks([1,10,100])
    xticklabels([1,10,100])
    ylim([1e-9,1e3])
    yticks(10.^[-9:6:3])
    gcaformat;

% Write to Source Data file.
% T = table(data.freq);
% T.Properties.VariableNames{1} = 'Frequency (Hz)';
% for i = 1:14
%     T{:,i+1} = round(log10(data.baseline(:,i)),3,'significant');
%     T.Properties.VariableNames{i+1} = ptIDs{i};
% end
% filename = 'E:\Research_Projects\004_Propofol\manuscript\Nature Communications\_final_submission\source_data.xlsx';
% writetable(T,filename,'Sheet','Figures 1i, 2b, 3f, 4g, S6a','Range','B2')