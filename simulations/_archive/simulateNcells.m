%%%%%%%%%%%% NYHM %%%%%%%%%%%%%%%%%
load('E:\Research_Projects\004_Propofol\Modelling\head_models\sa_nyhead.mat')
%%%%%%%%%%%% NYHM %%%%%%%%%%%%%%%%%

% neuron = 'NMO_68172'; 
neuron = 'L5E_oi15rpy4.swc';
% neuron = '86955';
propofol = 0; n = 10;
[err,prints] = system(['python "E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\code\simulation_alpha20220524.py" ' int2str(propofol) ' ' neuron ' ' int2str(n)]);
outputs = split(prints,char(10));
fileIdcs = find(cellfun(@(x)~isempty(strfind(x,'E:')),outputs));

for i = 1:length(fileIdcs)
    % rng(60);
    saveFile = outputs{fileIdcs(i)};
    data1(i) = load(saveFile);
    eeg1(:,i) = getEEG(data1(i).Q',sa);
    [p_baseline(:,i),f] = pmtm(eeg1(:,i),2,[],1e3*16);
end

neuron = 'NMO_68172'; 
% neuron = '86955';
propofol = 1; n = 50;
[err,prints] = system(['python "E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\code\runsimulation_alpha.py" ' int2str(propofol) ' ' neuron ' ' int2str(n)]);
outputs = split(prints,char(10));
fileIdcs = find(cellfun(@(x)~isempty(strfind(x,'E:')),outputs));

for i = 1:length(fileIdcs)
    % rng(60);
    saveFile = outputs{fileIdcs(i)};
    data2(i) = load(saveFile);
    eeg2(:,i) = getEEG(data2(i).Q',sa);
    [p_propofol(:,i),f] = pmtm(eeg2(:,i),2,[],1e3*16);
end


figureNB;
    plot(f,nanmedian(p_baseline,2),'linewidth',1,'color','k'); hold on;
    plot(f,nanmedian(p_propofol,2),'linewidth',1,'color','r'); hold on;
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,300]);
    gcaformat;
    xticks([1,10,100]);

actionchime