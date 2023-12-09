load(fullfile(dataFolder,'simulations','simulated_spectra_parameter_sensitivity'));

neuronTypes = {'L23E_oi24rpy1';'L23I_oi38lbc1';'L23I_oi38lbc1';'L4E_53rpy1';'L4E_j7_L4stellate';'L4E_j7_L4stellate';'L4I_oi26rbc1';'L4I_oi26rbc1';'L5E_oi15rpy4';'L5E_j4a';'L5I_oi15rbc1';'L5I_oi15rbc1';'L6E_51_2a_CNG';'L6E_oi15rpy4';'L6I_oi15rbc1';'L6I_oi15rbc1'};
nrnAbundance = [26.8,3.2,4.3,9.5,9.5,9.5,5.6,1.5,4.9,1.3,0.6,0.8,14,4.6,1.9,1.9];
mTypeCount = accumarray(findgroups(neuronTypes),nrnAbundance);
mTypeP = cumsum(mTypeCount)/sum(mTypeCount);
mType = @(p) interp1(mTypeP,1:11,p,'next','extrap');

% Define distributions of parameters
nPars = 8;
lambda_E = @(p) cdf('logn',p,log(0.5),1);
lambda_I = @(p) cdf('logn',p,log(2.5),1);
tauE = @(p) cdf('unif',p,1,3.5);
tauI = @(p) cdf('unif',p,5,20);
gE = @(p) cdf('unif',p,0.2e-3,2e-3);
gI = @(p) cdf('unif',p,0.2e-3,2e-3);
erev = @(p) cdf('unif',p,-75,-45);
gleak = @(p) cdf('unif',log10(p),-5,log10(0.005));
samplePars = @(P) [lambda_E(P(:,1)), lambda_I(P(:,2)), tauE(P(:,3)), ...
            tauI(P(:,4)), gE(P(:,5)), gI(P(:,6)), erev(P(:,7)), ...
            gleak(P(:,8)), mType(P(:,9))];
P = samplePars(pars')';

bins = linspace(0,1,50);
for i = 1:nPars
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
    E(j,:) = mean(oof_1_40(pars(9,:)==j,:));
    pr(j) = sum(pars(9,:)==j)/length(pars(9,:));
end
V1(9,:) = sum((E-mean(E)).^2.*pr);

V_total = var(oof_1_40);

p_names = {'\lambda_E','\lambda_I','\tau_E','\tau_I','g_E','g_I','E_{L}','g_{L}','m'};
figureNB;
subplot(1,2,1);
    [xSorted,I] = sort(V1(:,1),'descend');
    bar(xSorted/V_total(1))
    xticks(1:nPars);
    xticklabels(p_names(I));
    title('Slope Sensitivity')
subplot(1,2,2);
    [xSorted,I] = sort(V1(:,2),'descend');
    bar(xSorted/V_total(2))
    xticks(1:nPars);
    xticklabels(p_names(I));
    title('Amplitude Sensitivity')

bins = linspace(0,1,20);
V2_offset = zeros(nPars,nPars);
V2_beta = zeros(nPars,nPars);
for i = 1:nPars
    for j = i+1:nPars
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
V2_offset_full = V2_offset + V2_offset' + eye(nPars).*V1(1:nPars,1);
V2_beta_full = V2_beta + V2_beta' + eye(nPars).*V1(1:nPars,2);

figureNB(9.2,5);
axes('Position',[0.1,0.19,0.25,0.5]);
    imagesc(tril(V2_beta_full/V_total(2)))
    axis square;
    xticks(1:nPars);
    yticks(1:nPars);
    CB = colorbar;
    CB.Position = [0.1+0.25+0.02,0.19,0.03,0.5];
    CB.Label.String = 'Sensitivity index';
    xticklabels(p_names)
    yticklabels(p_names)
    colormap(clrsPT.sequential(1e3))
    xax = get(gca,'xaxis');
    xax.TickLabelRotation = 0;
    set(gca,'CLim',[0,0.1])
    gcaformat;
axes('Position',[0.585,0.19,0.25,0.5]);
    imagesc(tril(V2_offset_full/V_total(1)))
    axis square;
    xticks(1:nPars);
    yticks(1:nPars);
    xticklabels(p_names)
    yticklabels(p_names)
    colormap(clrsPT.sequential(1e3))
    CB = colorbar;
    CB.Position = [0.585+0.25+0.02,0.19,0.03,0.5];
    CB.Label.String = 'Sensitivity index';
    xax = get(gca,'xaxis');
    xax.TickLabelRotation = 0;
    set(gca,'CLim',[0,0.15])
    gcaformat;