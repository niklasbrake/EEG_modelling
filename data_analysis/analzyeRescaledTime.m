load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\analyzed\Cz_multitaper_mean.mat')
load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\raw\timeInformation.mat','timeInfo');
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

tRescaled = linspace(-1.5,0.5,50);
pRescaled = zeros(size(P,1),length(tRescaled),14);
for i = 1:14
    tr = -time/infusionTime(i);
    for k = 1:length(tRescaled)-1
        idx1 = interp1(tr,1:length(tr),tRescaled(k),'previous');
        idx2 = interp1(tr,1:length(tr),tRescaled(k+1),'next');
        pRescaled(:,k,i) = nanmean(P(:,idx1:idx2,i),2);
    end
end


oofPars = zeros(2,size(pRescaled,2),size(pRescaled,3));
for i = 1:14
    idcs = find(~isnan(pRescaled(1,:,i)));
    [temp,oofFun] = getFOOOF(ff,10.^pRescaled(:,idcs,i),false);
    oofPars(:,idcs,i) = temp';
end

pars_pt = zeros(13,size(pRescaled,2),size(pRescaled,3));
savePath = 'E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\rescaled_Sep2022';
for i = 1:14
    idcs = find(~isnan(pRescaled(1,:,i)));
    [temp,synFun] = synDetrend(ff,10.^pRescaled(:,idcs,i));
    pars_pt(:,idcs,i) = temp';
    pars = pars_pt(:,:,i);
    save(fullfile(savePath,['pt' int2str(i) '_rescaled_' date '.mat']),'pars');
end

for i = 1:14
    for j = 1:size(pRescaled,2)
        % OOF(:,j,i) = oofFun(ff,oofPars(:,j,i));
        ptBL(:,j,i) = synFun(ff,pars_pt(1:4,j,i));
    end 
end

P_biexp = pRescaled-ptBL;
P_oof = pRescaled-log10(OOF);

for i = 1:14
    P_basenorm(:,:,i) = pRescaled(:,:,i)-log10(pre(freq<150,i));
    P_biexp_basenorm(:,:,i) = P_biexp(:,:,i)-nanmedian(P_biexp(:,tfilt<-1,i),2);
    P_oof_basenorm(:,:,i) = P_oof(:,:,i)-nanmedian(P_oof(:,tfilt<-1,i),2);
end

