[sa,X] = network_simulation_beluga.getHeadModel;

passive = load('E:\Research_Projects\005_Aperiodic_EEG\data\simulations\active\LFPy_passive\simulation_data.mat')
active = load('E:\Research_Projects\005_Aperiodic_EEG\data\simulations\active\LFPy_active\simulation_data.mat')

idcs = find(passive.time>100);

passive.time = passive.time(idcs);
passive.V = passive.V(idcs,:);
passive.dipoles = passive.dipoles(idcs,:,:);

active.time = active.time(idcs);
active.V = active.V(idcs,:);
active.dipoles = active.dipoles(idcs,:,:);

location = 48108;

passive.eeg = network_simulation_beluga.getEEG(passive.dipoles,sa,location);
active.eeg = network_simulation_beluga.getEEG(active.dipoles,sa,location);

fl = 10.^linspace(log10(f(1)),log10(f(end)),100);
P_pas_log = [];
P_act_log = [];
for i = 1:size(P_passive,2)
    P_pas_log(:,i) = interp1(f,P_passive(2:end,i),fl,'linear','extrap');
    P_act_log(:,i) = interp1(f,P_active(2:end,i),fl,'linear','extrap');
end

dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';
load(fullfile(dataFolder,'cortical_column_Hagen','morphology_segmentations.mat'))
nrnSegs.L23E_oi24rpy1




red = clrsPT.qualitative_CM.red;
figureNB(12,8);
axes('Position',[0.03,0.08,0.36,0.77]);
    x = nrnSegs.L23E_oi24rpy1.x;
    y = nrnSegs.L23E_oi24rpy1.y;
    z = nrnSegs.L23E_oi24rpy1.z;
    line(x',z','color','k')
    set(gca,'DataAspectRatio',[1,1,1]);
    line([50,150],[-175,-175],'color','k','LineWidth',1.5)
    axis off;
axes('Position',[0.46,0.57,0.45,0.34]);
    plot(passive.time,active.V(:,2),'color',red);
    hold on;
    plot(passive.time,passive.V(:,2),'color','k')
    xlim([100,2000])
    line([100,100]+300,[-40,-20],'color','k','LineWidth',1.5);
    line([120,320]+300,[-40,-40],'color','k','LineWidth',1.5);
    ylim([-65,10])
    axis off;
    gcaformat;
axes('Position',[0.48,0.11,0.1,0.34]);
    boxplotNB(1,k/range(passive.time)*1e3,red,15);
    set(get(gca,'xaxis'),'visible','off');
    ylabel('Firing rate (Hz)');
    gcaformat;
axes('Position',[0.7,0.14,0.22,0.32]);
    plotwitherror(fl,P_pas_log,'M','color','k');
    hold on;
    plotwitherror(fl,P_act_log,'M','color',red);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylabel(['Unitary EEG spectrum (' char(956) 'V^2/Hz)'])
    xlabel('Frequency (Hz)');
    ylim([1e-25,1e-15])
    gcaformat;
    xlim([0.5,1e4])
    xticks([1,10,100,1000,1e4])
gcaformat(gcf);