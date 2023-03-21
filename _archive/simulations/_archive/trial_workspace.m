propofol = 0;
neuron = 'L5E_oi15rpy4.swc';
n=1;
xCorr = 0.95;

[err,prints] = system(['python "E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\code\simulationCorrInputs_update20220526.py" ' int2str(propofol) ' ' neuron ' ' int2str(n) ' ' num2str(xCorr,3)]);
returns = find(prints==char(10));
saveFile = prints(returns(end-1)+1:returns(end)-1)
data = load(saveFile);

fig = figureNB;
ax = subplot(1,3,1);
    for i = 1:length(data.morphology.segments)
        x = data.morphology.coordinates(data.morphology.segments{i},3:5);
        plot3(x(:,1),x(:,2),x(:,3),'k','parent',ax); hold on;
    end
    set(gca,'DataAspectRatio',[1,1,1])
    plot3(data.morphology.coordinates(:,3),data.morphology.coordinates(:,4),data.morphology.coordinates(:,5),'.')
    axis off;

ax = subplot(2,3,2);
    plot(data.t,data.V_soma,'color','r')
    gcaformat;
    xlabel('Time (ms)')
    ylabel('Voltage (mV)');


eeg = getEEG(data.Q',sa);

ax = subplot(2,3,5);
    plot(data.t,eeg,'color','k');
    gcaformat;
    xlabel('Time (ms)')
    ylabel('Voltage (uV)');

[P,f] = pmtm(eeg,2,[],1e3*16);; hold on;

ax = subplot(1,3,3);
    plot(f,P,'LineWidth',1,'color','k'); hold on;
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,300]);
    gcaformat;
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel(['Power (' char(956) 'V^2/Hz)'])