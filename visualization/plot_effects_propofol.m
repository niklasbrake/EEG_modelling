[sa,X] = network_simulation.getHeadModel;


folder = 'C:\Users\brake\Documents\temp\best';

F = dir(folder);
F = F(3:end);
for i=  1:length(F)
    data = load(fullfile(folder,F(i).name,'data.mat'));
    network(i) = data.network.importResults;;

    for j = 1:size(network(i).results.dipoles,3)
        dp(:,:,j) = resample(network(i).results.dipoles(2:end,:,j),1024,16e3);
    end
    qq{i} = gpuArray(dp);
end


tic
t = network(1).results.t;
idcs = sa.cortex2K.in_from_cortex75K(randi(2e3,30,1));
for j = 1:length(idcs)
    for i = 1:length(F)
        eeg = network_simulation.getEEG(qq{i},sa,idcs(j));
        temp = mypmtm(eeg,1024,2);
       if(j==1)
            P{i} = zeros(size(temp));
        end
        P{i} = P{i}+temp;
    end
end
toc

for i = 1:length(F)
    P{i} = P{i}/length(idcs)*pi/2;
end

f = 0.5:0.5:512;
clrs = clrsPT.sequential(length(P)+3);
clrs = clrs(3:end,:);
figureNB;
subplot(2,1,1);
    for i = 1:length(P)
        plot(f,mean(P{i},2),'LineWidth',2,'color',clrs(i,:))
        hold on;
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    gcaformat;
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlim([0.5,100])
    ylim([1e-18,1e-13])

subplot(2,1,2);

[full_model,synFun] = fittingmodel;
tRescaled = linspace(-1.5,0.5,200);
folder = 'E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\rescaled_manual\fitted';
for i = 1:14
load(fullfile(folder,['pt' int2str(i) '_rescaled_28-Sep-2022.mat']));
params(:,:,i) = pars;
end
ff = 0.5:0.5:150;
for i = 1:14
for j = 1:size(params,2)
ptBL(:,j,i) = synFun(ff,params(:,j,i));
end
end
N = 7;
m = 3;
clrs = clrsPT.sequential(N+m);
clrs = clrs(1+m:end,:);
idcs = find(and(tRescaled>-1,tRescaled<0));
k = floor(length(idcs)/N);
for i = 1:N
    BL = mean(mean(ptBL(:,idcs(k*(i-1)+1:k*i),:),2),3);
    plot(ff,10.^BL,'color',clrs(i,:),'LineWidth',1); hold on;
end
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([0.5,100])
xticks([1,10,100])
ylim([1e-2,1e3])
gcaformat
xlabel('Frequency (Hz)')
ylabel(['PSD (' char(956) 'V^2/Hz)'])
colormap(clrs)
CB = colorbar('location','south');
CB.TickLabels = {'Infusion','LOC'};
CB.Label.String = 'Rescaled time';
CB.TickDirection = 'out';
CB.TickLength = 0.05;
CB.Ticks = [0.5/N,1-0.5/N];