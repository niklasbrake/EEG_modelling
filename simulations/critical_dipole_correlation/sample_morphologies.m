fid = fopen(fname,'w');
for i = 1:100
    network = network_simulation_beluga(folder);
    % Initialize post network
    network = network.initialize_postsynaptic_network(2);
    M = {network.neurons(:).mType};
    fprintf(fid,'%s,%s\n',M{1},M{2});
end
fclose(fid)