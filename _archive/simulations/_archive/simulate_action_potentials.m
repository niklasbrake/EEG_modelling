N = 600;
tmax = 2e3;
masterPath = 'C:\Users\brake\Documents\temp\action_potentials2';
network = network_simulation(masterPath);
network = network.addActiveChannels();
network = network.initialize_postsynaptic_network(N);
network = network.initialize_presynaptic_network(0,tmax);
network.form_connections;
network = network.simulate;