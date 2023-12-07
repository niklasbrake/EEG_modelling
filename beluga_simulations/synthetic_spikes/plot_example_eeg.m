[sa,X] = network_simulation_beluga.getHeadModel;

% Import all data
folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\synthetic_spikes\mCombo1';
load(fullfile(folder,'simulation','simulation_data.mat'));

i = sa.cortex2K.in_from_cortex75K(263);
eeg = network_simulation_beluga.getEEG(dipoles,sa,i);


x = time(32e3:80e3);
y = eeg(32e3:80e3,:);
y(36e3:48e3,:) = nan;

figureNB(5.25,1);
axes('Position',[0,0,1,1]);
    plot(x,y(:,1),'LineWidth',0.5,'color','k');
    hold on;
    plot(x,y(:,2),'LineWidth',0.5,'color','r');
    line([4450,4450],0.25e-6+[-0.5,0.5]*1e-6,'color','k','LineWidth',2)
    line([4550,4800],0.25e-6+[-0.75,-0.75]*1e-6,'color','k','LineWidth',2)
    axis off