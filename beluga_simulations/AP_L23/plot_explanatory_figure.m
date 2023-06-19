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



passive.time = passive.time-passive.time(1);
red = clrsPT.qualitative_CM.red;
figureNB(12,8);
axes('Position',[0.03,0.08,0.36,0.77]);
    x = nrnSegs.L23E_oi24rpy1.x;
    y = nrnSegs.L23E_oi24rpy1.y;
    z = nrnSegs.L23E_oi24rpy1.z;
    line(x',z','color','k')
    set(gca,'DataAspectRatio',[1,1,1]);
    line([-50,-150],[-175,-175],'color','k','LineWidth',1.5)
    axis off;
axes('Position',[0.44,0.46,0.52,0.42]);
    plot(passive.time,active.dipoles(:,1,2),'color',[0.3,0.3,0]);
    line([-150,-50],[0,0],'color','k','LineWidth',1);
    text(-175,0,'Qx','FontSize',7,'HorizontalAlignment','right');
    hold on;

    plot(passive.time+200,active.dipoles(:,2,2)+25,'color',[0,0.3,0.3]);
    line([-150,-50]+200,[0,0]+25,'color','k','LineWidth',1);
    text(-175+200,25,'Qy','FontSize',7,'HorizontalAlignment','right');

    plot(passive.time+400,active.dipoles(:,3,2)+50,'color',[0.3,0,0.3]);
    line([-150,-50]+400,[0,0]+50,'color','k','LineWidth',1);
    text(-175+400,50,'Qz','FontSize',7,'HorizontalAlignment','right');

    line([2300,2500],[-20,-20],'color','k','Linewidth',1.5);
    line([2200,2200],[-20,0],'color','k','Linewidth',1.5);

    ylim([-20,100]);
    xlim([-200,2.5e3])


axes('Position',[0.44,0.1,0.52,0.42]);
    plot(passive.time,active.eeg(:,2),'color',red);
    hold on;
    plot(passive.time,passive.eeg(:,2),'color','k');
    line([500,700],[2,2]*1e-6,'color','k','LineWidth',1.5)
    line([400,400],[2,3]*1e-6,'color','k','LineWidth',1.5)
    axis off;

gcaformat(gcf);


load('E:\Research_Projects\005_Aperiodic_EEG\data\simulations\bAP_unitary_response\unitary_waveform')
AP_t = [-499:500]/16;
AP_waveform = interp1(AP_t,mean(unitary_waveform,2),AP_t(1):0.1:AP_t(end));
AP_t = AP_t(1):0.1:AP_t(end);

figureNB(12,4);
axes('Position',[0.09,0.21,0.08,0.7]);
    boxplotNB(1,k/range(passive.time)*1e3,red,15);
    set(get(gca,'xaxis'),'visible','off');
    ylabel('Firing rate (Hz)');
    gcaformat;
axes('Position',[0.25,0.21,0.08,0.7]);
    plot(AP_t,AP_waveform-median(AP_waveform),'color','k','LineWidth',1);
    ylabel(['EEG (' char(956) 'V)'])
    set(get(gca,'xaxis'),'visible','off');
    line([-5,5],[-2,-2]*1e-7,'color','k','LineWidth',2)
    text(0,-2.5e-7,'10 ms','FontSize',6,'HorizontalAlignment','center','VerticalAlignment','top');
    ylim([-2,15]*1e-7);
    xlim([-10,10]);
    gcaformat;
axes('Position',[0.44,0.22,0.228,0.7]);
    plotwitherror(f,P_passive(2:end,:),'CI');
    hold on;
    plotwitherror(f,P_active(2:end,:),'CI');
    plotwitherror(f0,psd0,'CI','color','k');
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlabel('Frequency (Hz)');
    xlim([0.5,8e3])
    ylim([1e-25,1e-15])
    xticks([1,10,100,1000,1e4])
    xlim([0.5,1e4])
    gcaformat;
axes('Position',[0.79,0.22,0.18,0.7]);
    x = mean(psd0+P_passive(2:end,:),2);
    y = mean(P_active(2:end,:),2);
    s = scatter(x(:),y(:),3,x(:)*[0,0,0],'filled');
    s.MarkerFaceAlpha = 0.2
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([1e-25,1e-15])
    ylim([1e-25,1e-15])
    yticks([1e-25,1e-15])
    xticks([1e-25,1e-15])
    line(get(gca,'xlim'),get(gca,'ylim'),'color','r','LineWidth',1)
    ylabel('EEG PSD (active)')
    xlabel('EEG PSD (passive \perp AP)')
    gcaformat;