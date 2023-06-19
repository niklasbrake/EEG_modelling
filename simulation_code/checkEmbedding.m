function [xSyn,mData,ei] = checkEmbedding(network)

    mType = network.neurons(1).mType;
    mTypeSegmentationData = fullfile(network_simulation_beluga.resourceFolder,'cortical_column_Hagen','morphology_segmentations.mat');
    load(mTypeSegmentationData)
    mData = nrnSegs.(mType);

    X = [mean(mData.x,2),mean(mData.y,2),mean(mData.z,2)];
    % cons = csvread('E:\Research_Projects\004_Propofol\data\simulations\raw\osc_2ndOrder\postsynaptic_network\connections.csv');
    cons = csvread(fullfile(network.postNetwork,'connections.csv'));
    idcs = find(cons(:,1)==1);

    xSyn = X(cons(idcs,2),:);
    C = cons(idcs,3);

    [ids,ts,ei] = network.getprenetwork(network.spikingFile);
    ei = ei(C);

    figureNB;
    line(mData.x',mData.y',mData.z','color','k');
    hold on;
    scatter3(xSyn(:,1),xSyn(:,2),xSyn(:,3),5,C,'filled');
    view([0,0]);
    set(gca,'DataAspectRatio',[1,1,1]);

end