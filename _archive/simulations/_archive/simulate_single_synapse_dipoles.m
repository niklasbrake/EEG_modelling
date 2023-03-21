M = 1100;
mTypes = randi(11,M,1);

masterPath = 'E:\Research_Projects\004_Propofol\data\simulations\raw\synapse_orientations';
network = network_simulation(masterPath);
network = network.initialize_postsynaptic_network(M,mTypes);
network.branchingIdx = 0;
EI = rand(1,M)>0.85;
network = network.setsynapsecount(M);
network.save_presynaptic_network((1:M)',150*ones(M,1),EI,network.getsynapsecount,fullfile(network.preNetwork,'spikeTimes.csv'));


load('E:\Research_Projects\004_Propofol\data\resources\cortical_column_Hagen\segment_areas.mat');

SD = zeros(M,3);
for j = 1:M
    mData = nrnSegs.(network.neurons(j).mType);
    sa = mData.area;
    seg = randperm(length(sa),1);
    pos = [mean(mData.x,2),mean(mData.y,2),mean(mData.z,2)];
    segDir = pos./vecnorm(pos,2,2);
    SD(j,:) = segDir(seg,:);
    filename = fullfile(network.postNetwork,'connections',[network.neurons(j).name '.csv']);
    fid = fopen(filename,'w');
    fprintf(fid,'%s,%s\n',int2str(seg),int2str(j));
    fclose(fid);
end
save(fullfile(masterPath,'synapse_compartments.mat'),'SD');

network.tmax = 240;
network = network.simulate;
