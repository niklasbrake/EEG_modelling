addpath('/home/nbrake/aperiodic_EEG_modelling/simulations/functions');
tmax = 2e3;

tau_decay =  [10,13,15,20,23];
tau_change = [2,3];
gam_change = [1,2,3];
lam_change = [1,0.9,0.8,0.7];

[tau,dtau,dgam,dlam] = ndgrid(tau_decay,tau_change,gam_change,lam_change);

param_changes = [tau(:),dtau(:),dgam(:),dlam(:)];
param_init = cat(2,tau_decay(:),ones(length(tau_decay),3));
params = [param_init;param_changes];

N = 200;

masterPath = '/home/nbrake/data/simulations/raw/parameter_search';
network = network_simulation_beluga(fullfile(masterPath,'_template'));
network = network.initialize_postsynaptic_network(N);
network = network.initialize_presynaptic_network(0,tmax);
network.form_connections;

for j = 54:size(params,1)
    name = ['param_combo_' int2str(j)];
    folder = fullfile(masterPath,name);

    network2 = network.copy_network(folder);

    fid = fopen(fullfile(folder,'params.csv'),'w');
    fprintf(fid,'tau,dtau,dgam,dlam\n');
    fprintf(fid,'%d,%d,%d,%.1f',params(j,:));
    fclose(fid);

    propofol = sum(params(j,1:3).*[100,10,1]);
    network2 = network2.addPropofol(propofol);

    rate = params(j,4);
    network2 = network2.setFiringRate(rate,rate/0.15);
    network2 = network2.initialize_presynaptic_network(0,tmax);
    network2.form_connections;
    system(['python functions/prep_simulations.py ' network2.postNetwork]);
    network2.save();
    % network2 = network2.simulate();
end