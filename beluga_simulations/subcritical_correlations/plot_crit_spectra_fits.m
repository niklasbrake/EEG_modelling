
mResults = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\analzye_simulations_crit.mat');
f = mResults.f;
clrs = clrsPT.sequential(10); clrs = clrs(5:end,:);


m = sort([0.99,0.98,0.95,0.86,0.63,0]);
tau0 = [11,13,15,22,100,100];
fig = figureNB(18.3,12);
for i = 1:6
    subplot(2,3,i);
    P = mean(mResults.P1(:,:,i),2);
    offset = P(1);
    P = P./offset;
    if(i>3)
        idcs = find(and(f>0,f<40));
        [params,synFun,full_model] = synDetrend(f(idcs),P(idcs),0,'lorenz',[tau0(i)*1e-3,0.01,0.02,-1]);
    else
        idcs = find(and(f>0,f<100));
        [params,synFun,full_model] = synDetrend(f(idcs),P(idcs),0,'lorenz',[tau0(i)*1e-3,0.003,0.02,-1]);
    end

    h = plotwitherror(f,mResults.P1(:,:,i),'Q','color',clrs(i,:),'LineWidth',1)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    hold on;
    plot(f,10.^full_model(f,params).*offset,'-k','LineWidth',1);

    plot(f,10.^full_model(f,params.*[1,1,1,0]+[0,0,0,-Inf]).*offset,'--k','LineWidth',0.5);
    plot(f,10.^full_model(f,params.*[1,1,0,1]+[0,0,-Inf,0]).*offset,'--k','LineWidth',0.5);
    xlim([1,40]);
    xticks([1,40]);

    tau(i) = params(1);
    title(['\tau_1  = ' int2str(1e3*tau(i)) ' ms'],'FontSize',8);
    ylim(10.^[-17,-14]);
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    L = legend(h,['m = ' num2str(m(i))]);
    L.Box = 'off';
    L.ItemTokenSize = [10,10];
    L.FontSize = 7
    
    gcaformat;
    drawnow;
end