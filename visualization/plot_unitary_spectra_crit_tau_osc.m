tauResults = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\analzye_simulations_tau.mat');
mResults = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\analzye_simulations_crit.mat');
f = mResults.f;

clrs = clrsPT.sequential(10); clrs = clrs(5:end,:);

figureNB(11.5,7);
m = sort([0.99,0.98,0.95,0.86,0.63,0]);
axes('Position',[0.11,0.65,0.15,0.32]);
    h=[];
    for i = 1:length(m)
        y = mResults.P1(:,:,i);
        % y = y./y(f==100,:);
        h(i) = plotwitherror(f,y,'CI','color',clrs(i,:));
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    % ylim(10.^[-18,-14])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    % str = sprintf('''m=%.2f'',',m(:));
    % Lstr = eval(['{' str(1:end-1) '}']);
    % L = legend(h,Lstr);
    L = legend(h,num2str(m(:)));
    L.Title.String = 'm value';
    L.Box = 'off';
    L.ItemTokenSize = [5,5];
    L.Position = [0.27,0.69,0.08,0.25];
    % text(5e3,1e-16,'Branching no. (m)','Rotation',-90','fontSize',7,'FontWeight','bold','HorizontalAlignment','center');
    gcaformat;

tau = [10,15,20,25,30,35];
axes('Position',[0.5,0.65,0.15,0.32]);
    h=[];
    for i = 1:length(tau)
        h(i) = plotwitherror(f,tauResults.P1(:,:,i),'CI','color',clrs(i,:));
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylim(10.^[-18,-14])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    % str = sprintf('''%d'',',tau(:));
    % Lstr = eval(['{' str(1:end-1) '}']);
    % L = legend(h,Lstr);
    L = legend(h,num2str(tau(:)));
    L.Box = 'off';
    L.Title.String = '\tau (ms)';
    L.ItemTokenSize = [5,5];
    L.Position = [0.66,0.69,0.08,0.25];
    % text(5e3,1e-16,'\tau_{decay} (ms)','Rotation',-90','fontSize',7,'FontWeight','bold','HorizontalAlignment','center');
    gcaformat;



tauResults = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\analzye_simulations_tau_osc.mat');
tauResults.P1 = tauResults.P1(:,:,[1,end]);
mResults = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\analzye_simulations_crit_osc.mat');
mResults.P1 = mResults.P1(:,:,[1,end]);
f = mResults.f;
% clrs = [clrsPT.qualitative_CM.blue;clrsPT.qualitative_CM.red];
clrs = clrs([1,end],:);

m = [0,0.99];
axes('Position',[0.11,0.11,0.15,0.32]);
    h=[];
    for i = 1:length(m)
        y = mResults.P1(:,:,i);
        % y = y./y(f==100,:);
        h(i) = plotwitherror(f,y,'CI','color',clrs(i,:));
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    % ylim(10.^[-18,-14])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    % str = sprintf('''m=%.2f'',',m(:));
    % Lstr = eval(['{' str(1:end-1) '}']);
    % L = legend(h,Lstr);
    L = legend(h,num2str(m(:)));
    L.Title.String = 'm value';
    L.Box = 'off';
    L.ItemTokenSize = [5,5];
    L.Position = [0.27,0.15,0.08,0.25];
    % text(5e3,1e-16,'Branching no. (m)','Rotation',-90','fontSize',7,'FontWeight','bold','HorizontalAlignment','center');
    gcaformat;

tau = [10,35];
axes('Position',[0.5,0.11,0.15,0.32]);
    h=[];
    for i = 1:length(tau)
        h(i) = plotwitherror(f,tauResults.P1(:,:,i),'CI','color',clrs(i,:));
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylim(10.^[-18,-14])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    % str = sprintf('''%d'',',tau(:));
    % Lstr = eval(['{' str(1:end-1) '}']);
    % L = legend(h,Lstr);
    L = legend(h,num2str(tau(:)));
    L.Box = 'off';
    L.Title.String = '\tau (ms)';
    L.ItemTokenSize = [5,5];
    L.Position = [0.66,0.15,0.08,0.25];
    % text(5e3,1e-16,'\tau_{decay} (ms)','Rotation',-90','fontSize',7,'FontWeight','bold','HorizontalAlignment','center');
    gcaformat;




P = squeeze(mean(tauResults.P1,2));
idcs = find(and(f>2,f<15));
P = P./P(1,:);

figureNB;

[params1,synFun,full_model] = synDetrend(f(idcs),P(idcs,1),0,'lorenz',[0.01,0.003,0,-1.2]);
subplot(2,2,1);
    plot(f,P(:,1))
    hold on;
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    plot(f,10.^synFun(f,params1))
    plot(f,10.^full_model(f,params1))
    xlim([2,100]);
subplot(2,2,3);
    plot(f,P(:,1)./10.^synFun(f,params1)')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([2,100]);

[params2,synFun,full_model] = synDetrend(f(idcs),P(idcs,2),0,'lorenz',[0.05,0.003,0.5,-1.6]);
subplot(2,2,2);
    plot(f,P(:,2))
    hold on;
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    plot(f,10.^synFun(f,params2))
    plot(f,10.^full_model(f,params2))
    xlim([2,100]);
subplot(2,2,4);
    plot(f,P(:,2)./10.^synFun(f,params2)')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([2,100]);