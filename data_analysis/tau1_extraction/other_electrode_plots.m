elName = {'Fz', 'Cz', 'Pz', 'C3', 'C4', 'CP3', 'CP4'};
saveFolder = 'E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\electrodes_rescaled';

figureNB(14,6);
for i = 1:7
    load(fullfile(saveFolder,sprintf('electrode%d_%s.mat',i,elName{i})));
    subplot(2,4,i);
    plotwitherror(linspace(-1.5,0.5,200),squeeze(pars(:,1,:)*1e3),'CI','color','k');
    xlim([-1.5,0.5]); ylim([10,75]);
    xlabel('Rescaled time');
    ylabel('\tau_1 (ms)')
    title(elName{i})
    xticks([-1,0]);
    xticklabels({'Infusion','LOC'});
    gcaformat;
end
gcaformat(gcf)