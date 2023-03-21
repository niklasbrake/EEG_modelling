clrs = clrsPT.lines(4);
blue = clrsPT.qualitative_CM.blue;

low = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\low_resolution_fits.mat');
high = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\high_resolution_fits.mat');
[coeff,score_low] = pca((low.X-mean(low.X))./std(low.X));
[~,score_high] = pca((high.X-mean(high.X))./std(high.X));

%%%%%%% Realisitic spectra
load('C:\Users\brake\Documents\temp\analzye_simulations_2.mat');
P_baseline = 1.2*P1(:,:,5);
asynch_low_res = 2*P_baseline;
osc_low_res = (P_baseline+P1(:,:,6));
crit_low_res = (P_baseline*1.25+0.75*P1(:,:,3));
tau_low_res = (P_baseline*1.25+0.75*P1(:,:,2));
f_low_res = f;

y1 = crit_low_res(f_low_res<50,:)./mean(asynch_low_res(1,:),2);
y2 = osc_low_res(f_low_res<50,:)./mean(asynch_low_res(1,:),2);
x = f_low_res(f_low_res<50);
y1 = mean(y1.*(1+0.3*randn(size(y1))),2);
y2 = mean(y2.*(1+0.3*randn(size(y2))),2);
% params3 = synDetrend(x,y1,1,'avalanches',[0.015    0.003   -0.37   -1.2    0.13    0    2    0.0000    0.2500]);
params3 = [0.012    0.0060   -0.05   -0.6    0.1    8    2.5    0.2    0.9];
params4 = [0.011    0.0060   -0.14   -0.4    0.15    4    3    0.5    0.7];
% params4 = synDetrend(x,y2,1,'avalanches',[0.015    0.003   -0.37   -1.2    0.13    0    2    0.0000    0.2500]);



load('C:\Users\brake\Documents\temp\analzye_simulations_10.mat')
f_high_res = f;
asynch_high_res = 2*mean(P(:,:,6),2);
osc_high_res = mean(P(:,:,6)+P(:,:,7),2); osc_high_res = osc_high_res./asynch_high_res(1);
crit_high_res = mean(P(:,:,6)*1.25+0.75*P(:,:,4),2); crit_high_res = crit_high_res./asynch_high_res(1);
asynch_high_res = asynch_high_res./asynch_high_res(1);

all_high_res = osc_high_res+crit_high_res;
all_high_res = all_high_res./all_high_res(1);
crit_high_res = crit_high_res./crit_high_res(1);
osc_high_res = osc_high_res./osc_high_res(1);

params1 = [0.015    0.003   -0.37   -1.2    0.13    9    2    0.0000    0.2500];
params2 = [0.0115    0.003   -0.04   -0.75    0.13    0    2    3    0.24];
[full_model,synFun,apFun] = fittingmodel('avalanches');

figureNB(13,8);

%% High resolution simulations %%%
subplot(2,3,1);
    plot(f_high_res,crit_high_res,'k','LineWidth',1);
    hold on;
    plot(f_high_res,10.^full_model(f_high_res,params1),'color',clrs(2,:),'LineWidth',1);
    plot(f_high_res,10.^apFun(f_high_res,params1),'color',clrs(3,:),'LineWidth',1);
    plot(f_high_res,10.^synFun(f_high_res,params1),'color',clrs(1,:),'LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.1,100])
    xticks(10.^[-1:2]);
    xticklabels([0.1,1,10,100]);
    ylim([0.01,10]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    title('+avalanches (m=0.98)','FontSize',7,'color',clrs(3,:),'FontWeight','normal');
    T = text(1.2,2.86,{'Avalanche','timescale'},'FontSize',6);
subplot(2,3,2);
    plot(f_high_res,osc_high_res,'k','LineWidth',1);
    hold on;
    plot(f_high_res,10.^full_model(f_high_res,params2),'color',clrs(2,:),'LineWidth',1);
    plot(f_high_res,10.^apFun(f_high_res,params2),'color',clrs(3,:),'LineWidth',1);
    plot(f_high_res,10.^synFun(f_high_res,params2),'color',clrs(1,:),'LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.1,100])
    xticks(10.^[-1:2]);
    xticklabels([0.1,1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    ylim([0.01,10]);
    yticks([]);
    title('+2 Hz oscillation','FontSize',7,'color',clrs(2,:),'FontWeight','normal');
    T = text(12.2,4,{'Peak'},'FontSize',6);
    T = text(0.77,0.11,{'Synaptic','timescale'},'FontSize',6);
subplot(2,3,3);
    for i = 1:4
        idcs = find(high.g==i);
        h(i) = scatter(score_high(idcs,1),score_high(idcs,2),30,clrs(i,:),'filled');
        hold on;
    end
    xlabel('PC1'); ylabel('PC2'); gcaformat;
    T1=text(-4,-2.1,{'Baseline','(m=0)'},'FontSize',6,'color',clrs(1,:));
    T2=text(-4.4,2.8,'+2 Hz oscillation','FontSize',6,'color',clrs(2,:));
    T3=text(-2.5,-4.6,'+avalanches (m=0.98)','FontSize',6,'color',clrs(3,:));
    T4=text(-0.5,4,'High tau (\tau=30 ms)','FontSize',6,'color',clrs(4,:));
    title('Freq. res. = 0.1 Hz','FontWeight','normal');
    ylim([-6,4.75])
    xlim([-5,5]);

subplot(2,3,4);
    plot(x,y1,'k','LineWidth',1);
    hold on;
    plot(x,10.^full_model(x,params3),'color',clrs(2,:),'LineWidth',1);
    plot(x,10.^apFun(x,params3),'color',clrs(3,:),'LineWidth',1);
    plot(x,10.^synFun(x,params3),'color',clrs(1,:),'LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.5,50])
    xticks([0.5,5,50]);
    xticklabels([0.5,5,50]);
    ylim([0.1,3.25])
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
subplot(2,3,5);
    plot(x,y2,'k','LineWidth',1);
    hold on;
    plot(x,10.^full_model(x,params4),'color',clrs(2,:),'LineWidth',1);
    plot(x,10.^apFun(x,params4),'color',clrs(3,:),'LineWidth',1);
    plot(x,10.^synFun(x,params4),'color',clrs(1,:),'LineWidth',1);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([0.5,50])
    xticks([0.5,5,50]);
    xticklabels([0.5,5,50]);
    xlabel('Frequency (Hz)');
    ylabel('log power');
    ylim([0.1,3.25])
    yticks([]);
subplot(2,3,6);
    for i = 1:4
        idcs = find(low.g==i);
        h(i) = scatter(score_low(idcs,1),score_low(idcs,2),30,clrs(i,:),'filled');
        hold on;
    end
    xlabel('PC1'); ylabel('PC2'); gcaformat;
    title('Freq. res. = 0.5 Hz','FontWeight','normal');
    ylim([-6,3]);
    xlim([-5,5]);
    % legend(h,{'Baseline (m=0)','+2 Hz oscillations','+avalanches (m=0.98)','\tau = 30 ms'});

gcaformat(gcf);



function computePCA
    %%%%%%% High resolution
    load('C:\Users\brake\Documents\temp\analzye_simulations_10.mat')
    f_high_res = f;
    asynch_high_res = 2*P(:,:,6);
    osc_high_res = P(:,:,6)+P(:,:,7);
    crit_high_res = P(:,:,6)*1.25+0.75*P(:,:,4);
    tau_high_res = 2*P(:,:,9);

    x = f_high_res(f_high_res<50);
    y0 = asynch_high_res(f_high_res<50,:)./asynch_high_res(100,:);
    y1 = osc_high_res(f_high_res<50,:)./osc_high_res(100,:);
    y2 = crit_high_res(f_high_res<50,:)./crit_high_res(100,:);
    y3 = tau_high_res(f_high_res<50,:)./tau_high_res(100,:);

    %%%%%%% Realisitic spectra
    load('C:\Users\brake\Documents\temp\analzye_simulations_2.mat');
    P_baseline = 1.2*P1(:,:,5);
    asynch_low_res = 2*P_baseline;
    osc_low_res = (P_baseline+P1(:,:,6));
    crit_low_res = (P_baseline*1.25+0.75*P1(:,:,3));
    tau_low_res = (P_baseline*1.25+0.75*P1(:,:,2));
    f_low_res = f;

    x = f_low_res(f_low_res<50);
    y0 = asynch_low_res(f_low_res<50,:)./mean(asynch_low_res(1,:),2);
    y1 = osc_low_res(f_low_res<50,:)./mean(asynch_low_res(1,:),2);
    y2 = crit_low_res(f_low_res<50,:)./mean(asynch_low_res(1,:),2);
    y3 = tau_low_res(f_low_res<50,:)./mean(asynch_low_res(1,:),2);
    y0 = y0.*(1+0.3*randn(size(y0)));
    y1 = y1.*(1+0.3*randn(size(y1)));
    y2 = y2.*(1+0.3*randn(size(y2)));
    y3 = y3.*(1+0.3*randn(size(y3)));

    for i = 1:14
        idcs = randi(14,14,1);
        y0_boot(:,i) = mean(y0(:,idcs),2);
        y1_boot(:,i) = mean(y1(:,idcs),2);
        y2_boot(:,i) = mean(y2(:,idcs),2);
        y3_boot(:,i) = mean(y3(:,idcs),2);
    end

    waitbar(0);
    P0 = [zeros(14,4),repmat([.14,0,2,0,0.25],14,1)];
    P1 = zeros(14,9);
    P2 = zeros(14,9);
    P3 = zeros(14,9);
    for i = 1:14
        waitbar(i/14);
        P0(i,1:4) = synDetrend(x,y0_boot(:,i),0,'lorenz',[0.012,0.002,0.05,-0.9]);
        P1(i,:) = synDetrend(x,y1_boot(:,i),1,'avalanches',P0(i,:));
        P2(i,:) = synDetrend(x,y2_boot(:,i),1,'avalanches',P0(i,:));
        P3(i,:) = synDetrend(x,y3_boot(:,i),1,'avalanches',[0.05,P0(i,2:end)]);
    end

    g = [ones(1,14),2*ones(1,14),3*ones(1,14),4*ones(1,14)];
    X = [P0;P1;P2;P3];
    X0 = X;
    [coeff,score] = pca((X0-mean(X0))./std(X0));

    clrs = clrsPT.lines(4);

    figureNB;
        scatter(score(:,1),score(:,2),30,clrs(g,:),'filled')
        xlabel('PC1'); ylabel('PC2'); gcaformat;

end