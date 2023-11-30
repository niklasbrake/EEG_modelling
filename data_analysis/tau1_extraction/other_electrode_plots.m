elName = {'Fz', 'Cz', 'Pz', 'C3', 'C4', 'CP3', 'CP4'}
figureNB(14,6);
for i = 1:7
    load(['E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\20230728_electrodes\rescaled_electrode' int2str(i) '.mat'])
    subplot(2,4,i);
    plotwitherror(rescaled.time,squeeze(rescaled.pars(1,:,:)*1e3),'CI','color','k');
    xlim([-1.5,0.5]); ylim([0,80]);
    xlabel('Rescaled time');
    ylabel('\tau_I (ms)')
    title(elName{i})
    xticks([-1,0]);
    xticklabels({'Infusion','LOC'});
    gcaformat;
end
gcaformat(gcf)