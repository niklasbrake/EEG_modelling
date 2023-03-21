load('E:\Research_Projects\004_Propofol\Modelling\head_models\sa_nyhead.mat')

neuron = '86955';

% figureNB;
for i = 1:9
    simulateEEG(0,neuron,sa)
    simulateEEG(1,neuron,sa)
end

actionchime

function simulateEEG(propofol,neuron,sa)

    filepath = 'E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\temp\';
    if(~propofol)
        filename = [filepath neuron '_simulation'];
    else
        filename = [filepath neuron '_propofol_simulation'];
    end

    system(['python "E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\code\runsimulation2.py" ' int2str(propofol) ' ' neuron])

    return;
    id = 1;
    file = [filename sprintf('%02d',id) '.mat'];
    while(exist(file))
        id = id+1;
        file = [filename sprintf('%02d',id) '.mat'];
    end
    id = id-1;
    file = [filename sprintf('%02d',id) '.mat'];

    data = load(file);

    subplot(2,3,1:2);
        plot(data.t,data.V(1,:)); hold on;
        gcaformat
    subplot(2,3,4:5);
        rng(10)
        eeg1 = getEEG(data.Q',sa);
        plot(data.t,eeg1); hold on;
        gcaformat

    subplot(1,3,3);
        p = [];
        rng(1)
        for i = 1:100
            eeg=getEEG(data.Q',sa);
            [p(:,i),f] = pmtm(eeg(data.t>500),2,[],1e3*16);; hold on;
        end
        p = nanmean(p,2);
        plot(f,p)
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        xlim([0.5,300]);
        gcaformat;
        xticks([1,10,100]);
    end