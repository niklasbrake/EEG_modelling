elName = {'Fz', 'Cz', 'Pz', 'C3', 'C4', 'CP3', 'CP4'};
saveFolder = 'E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\electrodes_rescaled';
filename = 'E:\Research_Projects\004_Propofol\manuscript\Nature Communications\_final_submission\source_data.xlsx';

% r = {'A','P','AE','AT','BI','BW','CK'};

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

    % Write Source Data files
    % tau = squeeze(pars(:,1,:)*1e3);
    % time = linspace(-0.5,1.5,200);
    % T = table(time(:));
    % T.Properties.VariableNames{1} = 'Time (rescaled)';
    % for k = 1:14
    %     T{:,k+(i==1)} = round(tau(:,k),3,'significant');
    %     T.Properties.VariableNames{k+(i==1)} = sprintf('pt. %d',k);
    % end
    % writetable(T,filename,'Sheet','Figure S7','Range',sprintf('%s2',r{i}))

end
gcaformat(gcf)