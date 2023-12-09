load(fullfile(dataFolder,'EEG_data','data_time_information.mat'));
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

rescaled = load(fullfile(dataFolder,'EEG_data','electrode2_Cz_rescaled_time.mat'));
rescaled.psd = log10(rescaled.psd);

aligned = load(fullfile(dataFolder,'EEG_data','electrode2_Cz.mat'));

load(fullfile(dataFolder,'EEG_data','Eq6_fits','electrode2_Cz.mat'),'pars');
rescaled.synPars(1:4,:,:) = permute(pars(:,1:4,:),[2,1,3]);
freq = rescaled.freq;

[full_model,synFun] = fittingmodel('eq6');
for i = 1:14
    for j = 1:size(rescaled.psd,2)
        synFit(:,j) = synFun(freq,rescaled.synPars(1:4,j,i));
    end
    rescaled.syn_detrended(:,:,i) = rescaled.psd(:,:,i)-synFit;

    rescaled.pre(:,:,i) = rescaled.psd(:,:,i)-nanmedian(rescaled.psd(:,rescaled.time<-1,i),2);
    rescaled.syn(:,:,i) = rescaled.syn_detrended(:,:,i)-nanmedian(rescaled.syn_detrended(:,rescaled.time<-1,i),2);

    synFitAligned = interp1(-t0(i)*rescaled.time',synFit',aligned.time,'cubic')';
    aligned.syn(:,:,i) = 10* (log10(aligned.psd(:,:,i)) - synFitAligned);
    P0(:,i) = nanmean(aligned.syn(:,aligned.time<t0(i),i),2);
end

deltaIdx = find(and(freq>=1,freq<4));
alphaIdx = find(and(freq>=8,freq<15));
betaIdx = find(and(freq>=15,freq<30));
aligned.alpha_syn = squeeze(nanmean(aligned.syn(alphaIdx,:,:)));
aligned.beta_syn = squeeze(nanmean(aligned.syn(betaIdx,:,:)));
aligned.delta_syn = squeeze(nanmean(aligned.syn(deltaIdx,:,:)));
aligned.alpha_pre = squeeze(nanmean(aligned.pre(alphaIdx,:,:)));
aligned.beta_pre = squeeze(nanmean(aligned.pre(betaIdx,:,:)));
aligned.delta_pre = squeeze(nanmean(aligned.pre(deltaIdx,:,:)));

for i = 1:14
    aligned.delta_syn(:,i) = aligned.delta_syn(:,i)-nanmean(aligned.delta_syn(aligned.time<=t0(i),i));
    aligned.alpha_syn(:,i) = aligned.alpha_syn(:,i)-nanmean(aligned.alpha_syn(aligned.time<=t0(i),i));
    aligned.beta_syn(:,i) = aligned.beta_syn(:,i)-nanmean(aligned.beta_syn(aligned.time<=t0(i),i));

    aligned.delta_syn(:,i) = fillgaps(aligned.delta_syn(:,i),5e3);
    aligned.alpha_syn(:,i) = fillgaps(aligned.alpha_syn(:,i),1024);
    aligned.beta_syn(:,i) = fillgaps(aligned.beta_syn(:,i),512);

    aligned.delta_pre(:,i) = fillgaps(aligned.delta_pre(:,i),5e3);
    aligned.alpha_pre(:,i) = fillgaps(aligned.alpha_pre(:,i),1024);
    aligned.beta_pre(:,i) = fillgaps(aligned.beta_pre(:,i),512);
end

% Not enough data points for computing baseline
aligned.delta_syn(:,3) = nan*aligned.delta_syn(:,3);

blue = clrsPT.qualitative_CM.blue;
clrs = clrsPT.lines(3);
clrs(2:3,:) = flip(clrs(2:3,:));

aligned.alpha_pre = smoothdata(aligned.alpha_pre,'movmedian',40);
aligned.beta_pre = smoothdata(aligned.beta_pre,'movmedian',40);
aligned.delta_pre = smoothdata(aligned.delta_pre,'movmedian',40);
aligned.alpha_syn = smoothdata(aligned.alpha_syn,'movmedian',40);
aligned.beta_syn = smoothdata(aligned.beta_syn,'movmedian',40);
aligned.delta_syn = smoothdata(aligned.delta_syn,'movmedian',40);

figureNB(7,6);
subplot(3,2,1);
    plotwitherror(aligned.time,aligned.alpha_pre,'CI','LineWidth',1,'color',clrs(1,:));
    xlim([-240,60]);
    ylim([-2.75,10])
    text(-220,10,'\alpha (8-15 Hz)','FontSize',6,'color',clrs(1,:),'VerticalAlignment','top');
    ylabel('Power (dB)')
    line([-300,60],[0,0],'color','k');
subplot(3,2,3);
    plotwitherror(aligned.time,aligned.beta_pre,'CI','LineWidth',1,'color',clrs(2,:));
    xlim([-240,60]);
    ylim([-2.5,7])
    text(-220,7,'\beta (15-30 Hz)','FontSize',6,'color',clrs(2,:),'VerticalAlignment','top');
    ylabel('Power (dB)')
    line([-300,60],[0,0],'color','k');
subplot(3,2,5);
    plotwitherror(aligned.time,aligned.delta_pre,'CI','LineWidth',1,'color',clrs(3,:));
    xlim([-240,60]);
    ylim([-3.5,13])
    text(-220,13,'\delta (1-4 Hz)','FontSize',6,'color',clrs(3,:),'VerticalAlignment','top');
    ylabel('Power (dB)')
    line([-300,60],[0,0],'color','k');
subplot(3,2,2);
    plotwitherror(aligned.time,aligned.alpha_syn,'CI','LineWidth',1,'color',clrs(1,:));
    xlim([-240,60]);
    ylabel('Power (dB)')
    line([-300,60],[0,0],'color','k');
    ylim([-2,7]);
subplot(3,2,4);
    plotwitherror(aligned.time,aligned.beta_syn,'CI','LineWidth',1,'color',clrs(2,:));
    xlim([-240,60]);
    ylabel('Power (dB)')
    line([-300,60],[0,0],'color','k');
    % ylim([-5,10])
subplot(3,2,6);
    plotwitherror(aligned.time,aligned.delta_syn,'CI','LineWidth',1,'color',clrs(3,:));
    xlim([-240,60]);
    ylabel('Power (dB)')
    line([-300,60],[0,0],'color','k');
    ylim([-4,7])
gcaformat(gcf)


% Write to Source Data files
% filename = 'E:\Research_Projects\004_Propofol\manuscript\Nature Communications\_final_submission\source_data.xlsx';
% time = linspace(-0.5,1.5,200);
% T = table(aligned.time(:));
% T.Properties.VariableNames{1} = 'Time (aligned)';
% writetable(T,filename,'Sheet','Figure S6d','Range','A2')

% T = table();
% for i = 1:14
%     T{:,i} = round(aligned.delta_pre(:,i),3,'significant');
%     T.Properties.VariableNames{i} = sprintf('pt. %d',i);
% end
% writetable(T,filename,'Sheet','Figure S6d','Range','C2')

% T = table();
% for i = 1:14
%     T{:,i} = round(aligned.alpha_pre(:,i),3,'significant');
%     T.Properties.VariableNames{i} = sprintf('pt. %d',i);
% end
% writetable(T,filename,'Sheet','Figure S6d','Range','R2')

% T = table();
% for i = 1:14
%     T{:,i} = round(aligned.beta_pre(:,i),3,'significant');
%     T.Properties.VariableNames{i} = sprintf('pt. %d',i);
% end
% writetable(T,filename,'Sheet','Figure S6d','Range','AG2')

% T = table();
% for i = 1:14
%     T{:,i} = round(aligned.delta_syn(:,i),3,'significant');
%     T.Properties.VariableNames{i} = sprintf('pt. %d',i);
% end
% writetable(T,filename,'Sheet','Figure S6d','Range','AV2')

% T = table();
% for i = 1:14
%     T{:,i} = round(aligned.alpha_syn(:,i),3,'significant');
%     T.Properties.VariableNames{i} = sprintf('pt. %d',i);
% end
% writetable(T,filename,'Sheet','Figure S6d','Range','BK2')

% T = table();
% for i = 1:14
%     T{:,i} = round(aligned.beta_syn(:,i),3,'significant');
%     T.Properties.VariableNames{i} = sprintf('pt. %d',i);
% end
% writetable(T,filename,'Sheet','Figure S6d','Range','BZ2')