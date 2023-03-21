% load('E:\Research_Projects\004_Propofol\Modelling\head_models\sa_nyhead.mat')


saveOutput1 = 'E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\asynchronous\baseline\passive\20220707\output.mat';

saveOutput2 = 'E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\critical\baseline\passive\20220707\output.mat';


[asynch.f,asynch.P,asynch.eeg,asynch.C_V,asynch.C_eeg,asynch.M] = analyzeSimulation(saveOutput1,sa);
[crit.f,crit.P,crit.eeg,crit.C_V,crit.C_eeg,crit.M] = analyzeSimulation(saveOutput2,sa);


figureNB;
plot(linspace(0,8e3,8193),mean(asynch.M,2)); hold on;
plot(linspace(0,8e3,8193),mean(crit.M,2)); hold on;

[asynch.PN,f] = pmtm(sum(asynch.eeg,2),2,[],1e3*16);
[crit.PN,f] = pmtm(sum(crit.eeg,2),2,[],1e3*16);

figureNB;
subplot(1,2,1);
    plot(asynch.f,sum(asynch.P,2)); hold on;
    plot(asynch.f,asynch.PN);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,150]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
subplot(1,2,2);
    plot(crit.f,sum(crit.P,2)); hold on;
    plot(crit.f,crit.PN);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,150]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])

function [f,P,eeg,C_V,C_eeg,M] = analyzeSimulation(saveOutput,sa)
    data = load(saveOutput);
    nrns = fieldnames(data);
    n = length(nrns);
    for i = 1:n
        Q(:,:,i) = data.(nrns{i}).sim.dipole;
        V(:,i) = data.(nrns{i}).sim.V_soma;
        preSyns{i} = data.(nrns{i}).syn.preSyns;
        temp = data.(nrns{i}).morph.x(data.(nrns{i}).syn.IDs,:);
        x{i} = temp./vecnorm(temp,2,2);
    end

    C_V = zeros(n,n);
    for ii = 1:n
        for jj = ii+1:n
            C_V(ii,jj) = corr(V(:,ii),V(:,jj));
        end
    end
    C_V = C_V+C_V';


    for ii = 1:size(Q,3)
        eeg(:,ii) = getEEG(Q(:,:,ii),sa,500);
        [P(:,i),f] = pmtm(eeg(:,ii),2,[],1e3*16);
    end

    C_eeg = zeros(n,n);
    k = 1;
    for ii = 1:n
        for jj = ii+1:n
            C_eeg(ii,jj) = corr(eeg(:,ii),eeg(:,jj));
            [M(:,k),f2] = mscohere(eeg(:,ii),eeg(:,jj),[],[],[],1e3*16);
            k = k+1;
        end
    end
    C_eeg = C_eeg+C_eeg';
    f2(1)
    f2(end)
    mean(diff(f2))
end
