load('E:\Research_Projects\004_Propofol\Modelling\head_models\sa_nyhead.mat')
X.vertices = sa.cortex75K.vc;
X.faces= sa.cortex75K.tri;
%%%%%%%%%%%% NYHM %%%%%%%%%%%%%%%%%

folder = 'E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\cortical_column\baseline\20220518'
F = dir(folder);
F = F(3:end);
eeg_baseline = [];
Q = zeros(32001,3,length(F));
for i = 1:length(F)
    data = load(fullfile(folder,F(i).name));
    % Q = Q+data.Q';
    Q(:,:,i) = data.Q';

end
Q = [Q(:,1,:),-Q(:,3,:),Q(:,2,:)];

idcs = sa.cortex2K.in_from_cortex75K;
N2 = length(idcs);

SA = 180000; %mm^2
rho3 = 108662; % neruons per mm^3
h = 2.082; % cortical height
rho2 = rho3/h; % neurons per mm^2

area = SA/N2; % Surface area per vertex
neurons = area*rho2

% Uncorrelated dipoles
q = zeros(size(Q,1),size(Q,2));
V = zeros(100,1);
m=2;
h = waitbar(0); tic;
for r = 1:100
    waitbar(r/100,h,['Estimated time remaining: ' char(duration(0,0,toc/r*(100-r)))]);
    density(r) = r/area;
    eeg = zeros(32001,1);
    q = r*sum(Q(:,:,randperm(100,r)),3);
    q2 = phaseran(q,N2);
    for i = 1:N2
        eeg = eeg+getEEG(q2(:,:,i),sa,idcs(i));
    end
    V(r) = var(eeg);
end

% Uncorrelated dipoles
q = zeros(size(Q,1),size(Q,2));
V2 = zeros(100,1);
h = waitbar(0); tic;
for r = 1:100
    waitbar(r/100,h,['Estimated time remaining: ' char(duration(0,0,toc/r*(100-r)))]);
    eeg = zeros(32001,1);
    q = Q(:,:,randi(100));
        for i = 1:N2
        eeg = eeg+getEEG(q2(:,:,i),sa,idcs(i));
    end
    V2(r) = var(eeg);
end


figureNB;
    plot(density,V); hold on;
    FT = fitlm(log(density),log(V));
    plot(1:rho2,exp(FT.predict(log(1:rho2)')))
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Density (neurons/mm^2)')
    ylabel(['EEG variance (' char(956) 'V^2)'])

tic;
eeg_baseline2 = zeros(32001,1);
var_corr = zeros(N,1);
k = 1;
for j = 1:N
    eeg_baseline2 = eeg_baseline2+getEEG(Q,sa,idx(j));
    var_corr(j) = var(eeg_baseline2);
end
toc


dataPSD = load('E:\Research_Projects\004_Propofol\Experiments\scalp_EEG\analyzed_data\psd_2s_1.9.mat')


n  = (1:N);
figureNB;
    plot(n,var_corr); hold on;
    x = log(n); y = log(var_corr);
    FT = fitlm(x,y);
    plot(n,exp(FT.predict(log(n)')));
    set(gca,'xscale','log')
    set(gca,'yscale','log')

    fill([10e9,80e9,80e9,10e9],[1e-20,1e-20,1e2,1e2],'k','FaceAlpha',0.2,'LineStyle','none');
    fill([1000,80e9,80e9,1000],[1e-4,1e-4,1e2,1e2],'k','FaceAlpha',0.2,'LineStyle','none');
    xlim([1000,1e11]);
    ylim([1e-20,1e2]);

neurons = 25e9/100;

[p_corr,f] = pmtm(eeg_baseline2,2,[],1e3*16);
[p_rand,f] = pmtm(eeg_baseline_random,2,[],1e3*16);


for i = 1:1000
    [p_single(:,i),f] = pmtm(getEEG(Q,sa),2,[],1e3*16);
end
p_single = nanmedian(p_single,2);

fig = figure('color','w');
subplot(1,2,1);
    plot_mesh_brain(X); 
    hold on;
    plot3(X.vertices(idx,1),X.vertices(idx,2),X.vertices(idx,3),'.','color',[0.5,0.5,1],'MarkerSize',2)
subplot(1,2,2);
    % plot(f,p_single,'linewidth',1,'color','b');hold on;
    plot(f,neurons*p_single,'linewidth',1,'color','r');hold on;
    plot(f,neurons.^(1.5)*p_single,'linewidth',1,'color','r');hold on;
    plot(dataPSD.data.freq,nanmean(nanmean(dataPSD.data.psd(:,1000:2000,:),2),3),'LineWidth',1,'color','k')
    % plot(f,p_corr,'linewidth',1,'color','k'); hold on;
    % plot(f,p_rand,'linewidth',1,'color','r'); hold on;
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,300]);
    gcaformat;
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel(['Power (' char(956) 'V^2/Hz)'])