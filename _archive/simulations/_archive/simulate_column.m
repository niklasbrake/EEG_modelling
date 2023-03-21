load('E:\Research_Projects\004_Propofol\Modelling\head_models\sa_nyhead.mat')

N = 100; % Total size of simulation
%  Hagen et al., 2016
neurons = {'L23E_oi24rpy1.swc';'L23I_oi38lbc1.swc';'L23I_oi38lbc1.swc';'L4E_53rpy1.swc'; ...
            'L4E_j7_L4stellate.swc';'L4E_j7_L4stellate.swc';'L4I_oi26rbc1.swc';'L4I_oi26rbc1.swc';
            'L5E_oi15rpy4.swc';'L5E_j4a.swc';'L5I_oi15rbc1.swc';'L5I_oi15rbc1.swc';'L6E_51-2a.CNG.swc'; ...
            'L6E_oi15rpy4.swc';'L6I_oi15rbc1.swc';'L6I_oi15rbc1.swc'};

relOccurance = [26.8,3.2,4.3,9.5,9.5,9.5,5.6,1.5,4.9,1.3,0.6,0.8,14,4.6,1.9,1.9];
occuranceCDF = cumsum(relOccurance)/sum(relOccurance);
randM = interp1(occuranceCDF,1:length(occuranceCDF),rand(100,1),'next','extrap');
M = histcounts(randM,0.5:1:length(neurons)+0.5);

neurons = unique(neurons);
propofol = 0;
for i = 1:length(neurons)
    neuron = neurons{i};
    % n = M(i);
    n=2;
    [err,prints] = system(['python "E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\code\simulation_update20220526.py" ' int2str(propofol) ' ' neuron ' ' int2str(n)]);
end