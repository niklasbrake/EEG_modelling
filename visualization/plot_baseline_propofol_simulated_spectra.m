
data1 = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\unitary_spectrum_asynch.mat');

% data2 = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\unitary_spectrum_propofol_highAmplitude3.mat');
data2 = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\unitary_spectrum_propofol_doubleDouble.mat');



figureNB;
    plot(data1.f,mean(data1.P,2),'LineWidth',1)
    hold on;
    plot(data2.f,mean(data2.P,2),'LineWidth',1)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    gcaformat;
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlim([0.5,100])