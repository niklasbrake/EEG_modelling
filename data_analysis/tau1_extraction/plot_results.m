load(fullfile(dataFolder,'EEG_data','data_time_information.mat'));
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

load(fullfile(dataFolder,'EEG_data','electrode2_Cz.mat'),'psd','freq','time');

% Compute Baseline spectrum
for i = 1:14
    preInfusion(:,i) = nanmedian(psd(:,and(time>=t0(i)-10,time<t0(i)),i),2);
    preLOC(:,i) = nanmedian(psd(:,and(time>=-10,time<0),i),2);
end

iNoise = find(and(freq>55,freq<65));
preInfusion(iNoise,:) = nan; preInfusion = fillgaps(preInfusion,5);
preLOC(iNoise,:) = nan; preLOC = fillgaps(preLOC,5);

% Convert to Baseline-normalized decibels
dB = log10(nanmedian(psd,3));

red = clrsPT.qualitative_CM.red;
clrs = clrsPT.sequential(7);
blue = clrsPT.qualitative_CM.blue;

CM = clrsPT.iridescent(1e3);
black = CM(end,:);

fig = figureNB(8.9,10);
axes('Position',[0.067,0.76,0.91,0.2]);

axes('Position',[0.11,0.465,0.375,0.26]);
    imagesc(time,freq,dB);
    axis xy
    ylim([0.5,40])
    CB = colorbar;
    CB.Location = 'eastoutside';
    CB.Position(1) = 0.5;
    CB.Position(3) = 0.03;
    set(gca,'CLim',[-2,3])
    CM = clrsPT.iridescent(50);
    idcs = interp1(linspace(0,1,50),1:50,linspace(-0.3,1.5,50),'nearest','extrap');
    colormap(flip(CM(idcs,:)));
    set(gca,'fontsize',7);
    xlabel('LOC-aligned time (s)');
    ylabel('Frequency (Hz)');
    line([0,0],get(gca,'ylim'),'color','k','linestyle','-','linewidth',1);
    line([0,0],get(gca,'ylim'),'color','w','linestyle','--','linewidth',1);
    xlim([-300,60]);
    CB.Ticks = [-2,0,2,4];
    CB.TickLabels = {'10^{-2}','10^0','10^2','10^{4}'};
    CB.FontSize = 7;

axes('Position',[0.71,0.465,0.25,0.26]);
    plotwitherror(freq,log10(preInfusion),'CI','linewidth',1,'color','k');
    plotwitherror(freq,log10(preLOC),'CI','linewidth',1,'color',red);
    xlim([0.5,300])
    set(gca,'xscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlim([0.5,100]);
    gcaformat;
    xlabel('Frequency (Hz)');
    yl = ylabel(['PSD (' char(956) 'V^2/Hz)']);
    yl.Position(1:2) = [0.1,0.36];
    ylim([-2,3])
    yticks([-2,0,2,4]);
    yticklabels({'10^{-2}','10^0','10^2','10^{4}'});
    text(0.7,-0.5,'Baseline','fontsize',7);
    text(2.3,2.6,'Pre-LOC (0-10 s)','fontsize',7,'color',red);


load(fullfile(dataFolder,'EEG_data','Eq6_fits','electrode2_Cz_baseline_and_preLOC.mat'),'pars_preInfusion','pars_preLOC');

[full_model,synFun] = fittingmodel;
for i = 1:14
    pSynPre(:,i) = 10.^synFun(freq,pars_preInfusion(:,i));
    pSynPost(:,i) = 10.^synFun(freq,pars_preLOC(:,i));
end

ptIdx = 1;

synPre = pars_preInfusion(:,ptIdx);
synPost = pars_preLOC(:,ptIdx);

pars_preInfusion(1,ptIdx)*1e3
pars_preLOC(1,ptIdx)*1e3

axes('Position',[0.11,0.09,0.22,0.26]);
    plot(freq,preInfusion(:,ptIdx),'color','k','LineWidth',1);
    hold on;
    plot(freq,10.^synFun(freq,synPre),'color',blue,'linestyle',':','LineWidth',1);
    plot(freq,10.^full_model(freq,synPre),'linestyle','-','linewidth',1,'color',blue);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlim([0.5,100]);
    gcaformat;
    xl = xlabel('Frequency (Hz)');
    xl.Position(1:2) = [150,1.1e-4];
    yl = ylabel(['PSD (' char(956) 'V^2/Hz)']);
    yl.Position(1) = 0.1;
    ylim([1e-3,1e3])
    set(get(gca,'yaxis'),'MinorTick','off');
    text(0.6,3e-3,'Baseline (pt.1)','FontSize',7)
axes('Position',[0.36,0.09,0.22,0.26]);
    plot(freq,preLOC(:,ptIdx),'color','k','LineWidth',1);
    hold on;
    plot(freq,10.^synFun(freq,synPost),'color',blue,'linestyle',':','LineWidth',1);
    plot(freq,10.^full_model(freq,synPost),'linestyle','-','linewidth',1,'color',blue);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlim([0.5,100]);
    gcaformat;
    ylim([1e-3,1e3])
    set(get(gca,'yaxis'),'MinorTick','off');
    yticklabels({});
    text(0.6,3e-3,'Pre-LOC (pt.1)','FontSize',7,'color','k')


fits = load(fullfile(dataFolder,'EEG_data','Eq6_fits','electrode2_Cz.mat'));
fits.time = linspace(-1.5,0.5,200);
tau = squeeze(fits.pars(:,1,:))*1e3;

axes('Position',[0.71,0.09,0.25,0.26]);
    plotwitherror(fits.time,tau,'CI','color','k')
    xticks([-1,0]);
    xticklabels({'Infusion','LOC'})
    xlabel('Rescaled time')
    ylabel('\tau_1 (ms)')
    xlim([-1.4,0.4])
    ylim([10,65]);
    gcaformat


gcaformat(gcf)

labelpanel(0.01,0.95,'a',true);
labelpanel(0.01,0.7,'b',true);
labelpanel(0.61,0.7,'c',true);
labelpanel(0.01,0.325,'d',true);
labelpanel(0.61,0.325,'e',true);;

% Write to Source Data file
% filename = 'E:\Research_Projects\004_Propofol\manuscript\Nature Communications\_final_submission\source_data.xlsx';

% T = table(eeg_example.time(:),eeg_example.timedomain(:));
% T.Properties.VariableNames = {'Time (s)','EEG (μV)'};
% writetable(T,filename,'Sheet','Figure 8a','Range','B2')

% T = table(freq(freq<=40));
% T.Properties.VariableNames{1} = 'Frequency (Hz)';
% for i = 1:size(dB,2)
%     T{:,i+1} = round(dB(freq<=40,i),3,'significant');
%     T.Properties.VariableNames{i+1} = sprintf('%.1f s',time(i));
% end
% writetable(T,filename,'Sheet','Figure 8b','Range','B2')

% T = table(freq);
% T.Properties.VariableNames{1} = 'Frequency (Hz)';
% for i = 1:14
%     T{:,i+1} = round(log10(preInfusion(:,i)),3,'significant');
%     T.Properties.VariableNames{i+1} = sprintf('pt. %d',i);;
% end
% writetable(T,filename,'Sheet','Figure 8c, d','Range','B2')

% T = table(freq);
% T.Properties.VariableNames{1} = 'Frequency (Hz)';
% for i = 1:14
%     T{:,i+1} = round(log10(preLOC(:,i)),3,'significant');
%     T.Properties.VariableNames{i+1} = sprintf('pt. %d',i);;
% end
% writetable(T,filename,'Sheet','Figure 8c, d, S6b','Range','B2')