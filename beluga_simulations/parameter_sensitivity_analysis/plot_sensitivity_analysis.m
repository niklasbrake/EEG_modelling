
load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis\fitted_spectra.mat')

bins = linspace(0,1,50);
for i = 1:6
    binnedPar = interp1(bins,1:50,P(i,:),'nearest');
    E = zeros(50,2);
    pr = zeros(50,1);
    for j = 1:50
        E(j,:) = mean(oof_1_40(binnedPar==j,:));
        pr(j) = sum(binnedPar==j)/length(binnedPar);
    end
    V1(i,:) = sum((E-mean(E)).^2.*pr);
end

E = zeros(11,2);
pr = zeros(11,1);
for j = 1:11
    E(j,:) = mean(oof_1_40(pars(7,:)==j,:));
    pr(j) = sum(pars(7,:)==j)/length(pars(7,:));
end
V1(7,:) = sum((E-mean(E)).^2.*pr);

V_total = var(oof_1_40);

p_names = {'\lambda_E','\lambda_I','\tau_E','\tau_I','E_{L}','g_{L}','m'};
figureNB;
subplot(1,2,1);
    [xSorted,I] = sort(V1(:,1),'descend');
    bar(xSorted/V_total(1))
    xticks(1:7);
    xticklabels(p_names(I));
subplot(1,2,2);
    [xSorted,I] = sort(V1(:,2),'descend');
    bar(xSorted/V_total(2))
    xticks(1:7);
    xticklabels(p_names(I));



bins = linspace(0,1,20);
V2_offset = zeros(6,6);
V2_beta = zeros(6,6);
for i = 1:6
    for j = i+1:6
        binnedPar1 = interp1(bins,1:length(bins),P(i,:),'nearest');
        binnedPar2 = interp1(bins,1:length(bins),P(j,:),'nearest');
        idx = (binnedPar1-1)*length(bins)+binnedPar2;
        uIdx = ((1:20)-1)*20+(1:20)';
        uIdx = uIdx(:);
        E = zeros(length(uIdx),2);
        pr = zeros(length(uIdx),1);
        for k = 1:length(uIdx)
            E(k,:) = mean(oof_1_40(idx==uIdx(k),:));
            pr(k) = sum(idx==uIdx(k))/length(idx);
        end
        temp = sum((E-mean(E)).^2.*pr) - V1(i,:) - V1(j,:);
        V2_offset(i,j) = temp(1);
        V2_beta(i,j) = temp(2);
    end
end
V2_offset_full = V2_offset + V2_offset' + eye(6).*V1(1:6,1);
V2_beta_full = V2_beta + V2_beta' + eye(6).*V1(1:6,2);

figureNB(7,3.2);
axes('Position',[0.1,0.19,0.25,0.5]);
    imagesc(tril(V2_beta_full/V_total(2)))
    axis square;
    xticks(1:6);
    yticks(1:6);
    CB = colorbar;
    CB.Position = [0.1+0.25+0.02,0.19,0.03,0.5];
    xticklabels(p_names)
    yticklabels(p_names)
    colormap(clrsPT.sequential(1e3))
    xax = get(gca,'xaxis');
    xax.TickLabelRotation = 0;
    set(gca,'CLim',[0,0.1])
axes('Position',[0.585,0.19,0.25,0.5]);
    imagesc(tril(V2_offset_full/V_total(1)))
    axis square;
    xticks(1:6);
    yticks(1:6);
    xticklabels(p_names)
    yticklabels(p_names)
    colormap(clrsPT.sequential(1e3))
    CB = colorbar;
    CB.Position = [0.585+0.25+0.02,0.19,0.03,0.5];
    xax = get(gca,'xaxis');
    xax.TickLabelRotation = 0;
    set(gca,'CLim',[0,0.15])
gcaformat(gcf)