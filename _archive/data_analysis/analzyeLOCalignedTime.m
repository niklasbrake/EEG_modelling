load('E:\Research_Projects\004_Propofol\Experiments\scalp_EEG\analyzed_data\Cz_multitaper_mean.mat')
load('E:\Research_Projects\004_Propofol\Experiments\scalp_EEG\raw_data\timeInformation.mat','timeInfo');
infusionTime = timeInfo.infusion_onset-timeInfo.object_drop;
for i = 1:size(psd,3)
    pre(:,i) = nanmedian(psd(:,time<infusionTime(i),i),2);
    post(:,i) = nanmedian(psd(:,and(time>=0,time<60),i),2);
    for j = 1:floor(freq(end)/60)
        idcs = find(and(freq>60*j-5,freq<60*j+5));
        pre(idcs,i) = nan;
        post(idcs,i) = nan;
    end
    pre(:,i) = fillgaps(pre(:,i),5);
    post(:,i) = fillgaps(post(:,i),5);
end

psd(psd==0) = nan;

ff = freq(freq<150);

P = nan(299,size(psd,2),size(psd,3));
for idx = 1:14
    idx
    goodIdcs = find(~isnan(psd(1,:,idx)));
    t0 = goodIdcs(1); t1 = goodIdcs(end);
    y = log10(psd(freq<150,t0:t1,idx));
    goodIdcs = find(~isnan(y(1,:)));
    for i = 1:length(goodIdcs)
        for j = 1:2
            idcs = find(and(freq>60*j-5,freq<60*j+5));
            y(idcs,goodIdcs(i)) = nan;
        end
        y(:,goodIdcs(i)) = fillgaps(y(:,goodIdcs(i)),5);
    end
    for i = 1:size(y,1)
        y(i,:) = fillgaps(y(i,:),5);
    end
    P(:,t0:t1,idx) = y;
end

preP = permute(log10(pre(freq<150,:)),[1,3,2]);
figureNB;
    imagesc(time,ff,nanmedian(P,3));
    ylim([0,50])
    colormap('jet')
    % set(gca,'CLim',[-10,18])
    colorbar;
    axis xy
    xlim([-150,60])
    ylabel('Frequency (Hz)')
    xlabel('Time to LOC (s)')
    gcaformat;

idcs = find(and(time>-300,time<60));
ds = 20;
pfilt = zeros(size(P,1),ceil(length(idcs)/ds),size(P,3));
for i = 1:size(P,3)
    % for j = 1:size(P,1)
    %     pfilt(j,:,i) = filtfilt(b,a,P(j,idcs,i));
    % end
    pfilt(:,:,i) = resample(P(:,idcs,i)',1,ds)';
end
tfilt = downsample(time(idcs),ds);


oofPars = zeros(2,size(pfilt,2),size(pfilt,3));
for j = 1:14
    [temp,fun] = getFOOOF(ff,10.^pfilt(:,:,j),false);
    oofPars(:,:,j) = temp';
end

OOF = [];
for i = 1:size(A,1)
    for j = 1:size(A,2)
        OOF(:,i,j) = fun(ff,[A(i,j),B(i,j)]);
    end
end

savePath = 'E:\Research_Projects\004_Propofol\Modelling\data_fitting\data\LOC_time';
pars_pt = zeros(13,size(pfilt,2),size(pfilt,3));
for i = 1:14
    [temp,synFun] = synDetrend(ff,10.^pfilt(:,:,i));
    pars_pt(:,:,i) = temp';
    pars = pars_pt(:,:,i);
    save(fullfile(savePath,['pt' int2str(i) '_LOCaligned.mat']),'pars');
end

for i = 1:14
    fileName = ['pt' int2str(i) '_spectrogram'];
    myPath = pwd;
    file0 = strrep(fullfile(myPath,[fileName '.csv']),'\','/');
    csvwrite(file0,[ff,10.^pfilt(:,:,i)]);
    pyFunction = 'C:\Users\brake\Documents\GitHub\Propofol2022\functions\model_fitting\detrending.py';
    cmd = ['python ' pyFunction ' ' file0];
    [err,prints] = system(cmd);
    file{i} = file0;
end

for i = 1:14
    pars_pt(:,:,i) = csvread([file{i}(1:end-4) '_params.csv']);
end

[F,FBL] = fittingmodel;
for i = 1:14
    for j = 1:size(pfilt,2)
        y = pfilt(:,j,i);
        offset = y(1);
        y = y-offset;
        ptBL(:,j,i) = FBL(pars_pt(j,1:4,i),ff)+offset;
    end 
end

PP = pfilt-ptBL;


figureNB;
i = 2;
h(1) = plot(ff,pfilt(:,i,1)-pfilt(1,i,1)); hold on;
h(2) = plot(ff,FBL(pars_pt1(i,:),ff));
set(gca,'xscale','log');
ylim([-7,4]);
xlim([0.5,150]);
for i = 1:size(pars_pt1,1)
    h(1).YData=pfilt(:,i,1);
    pt1_BL(:,i) = FBL(pars_pt1(i,:),ff)+pfilt(1,i,1);
    h(2).YData=pt1_BL(:,i);
    pause(0.1);
    drawnow;
end



figureNB;
    plotwitherror(freq,preDetrended,'CI','color',[0.5,0.5,1]); hold on;
    plotwitherror(freq,postDetrended,'CI','color',[1,0.5,0.5]);
    set(gca,'xscale','log');
    xlim([0.5,150]);
    gcaformat;