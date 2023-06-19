[sa,X] = network_simulation_beluga.getHeadModel;

locations = randi(size(sa.cortex75K.vc,1),100,1); % Random location
load('E:\Research_Projects\004_Propofol\data\simulations\raw\aperiodic_sensitivity\analyzed_results.mat')

location = 64030;
eeg = network_simulation_beluga.getEEG(dipoles,sa,location);
eeg1 = resample(eeg,1e3,16e3);
[psd,f2] = pmtm(detrend(eeg1,'constant'),2,[],1e3);


load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\oscillation_interactions.mat');
f = 0.1:0.1:500;


figureNB(30,17);
for i = 1:4
    for j = 1:4
        psd0 = mean(P2(:,:,(i-1)*4+1),2);
        psd1 = mean(P2(:,:,(i-1)*4+j),2);
        subplot(5,4,4*(i-1)+j);
        if(i>1)
            subplot(5,4,4*(i)+j);
        end
            if(j==1)
                clrs = [0.6,0.6,0.6;0,0,0];
                psd0 = mean(P2(:,:,1),2);
            else
                clrs = [0,0,0;1,0,0];
            end
            plot(f,psd0,'color',clrs(1,:),'LineWidth',1);
            hold on;
            plot(f,psd1,'color',clrs(2,:),'LineWidth',1);
            set(gca,'xscale','log')
            set(gca,'yscale','log')
            xlim([0.5,100])
            xticks([1,100]);
            ylim([1e-17,1e-15]);
            if(i<3)
                xticklabels({});
            else
                xticklabels([1,100]);
            end
            % yticks([]);
            gcaformat;
            drawnow;
            xticks([]);
            yticks([]);
            box on;
            if(i==4);
                xticks([1,10,100]);
                xticklabels([1,10,100]);
                xlabel('Frequency (Hz)');
            end
    end
end


idcs(1,:) = (data.parameters(1,:) == 10);
idcs(2,:) = (data.parameters(2,:) == -58.5);
idcs(3,:) = (data.parameters(3,:) == 0.5);
idcs(4,:) = (data.parameters(4,:) == 0.98);

for i = 1:4
    if(i==4)
        subplot(5,4,5);
        clrs = [0.6,0.6,0.6;0,0,0];
    else
        subplot(5,4,5+i);
        clrs = [0,0,0;1,0,0];
    end
    K = unique(data.parameters(i,:));
    cla;
    for k = 1:2
        temp = setdiff(1:4,i);
        idcs0 = find(prod(idcs(temp,:)));
        idcs2 = find(data.parameters(i,:)==K(k));
        P_crit(:,i) = mean(psd(:,intersect(idcs0,idcs2)),2);
        plot(f2,P_crit(:,i),'color',clrs(k,:),'LineWidth',1);
        hold on;
    end
    gcaformat;
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    ylim([1e-17,1e-13]);
    xticks([]);
    yticks([]);
    box on;
end